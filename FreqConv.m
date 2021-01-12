classdef FreqConv < audioPlugin
    properties
        mix = 50;
        max_frames = 150;
    end
    properties (Access = private)
        input_buffer;
        output_buffer;
        overlap_buffer;
        delay_fft_buffer_left;
        delay_fft_buffer_right;
        IR_frames_left;
        IR_frames_right;
        n_IR_frames = 0;
        delay_buffer_size = 1000000;
        channels = 2;
        buffer_size = 512;
        write_head = 1;
        NFFT = 1024;
    end
    properties (Constant)
        PluginInterface = audioPluginInterface( ...
            audioPluginParameter("mix",...
            "Mapping",{"lin",0,100}),...
            audioPluginParameter("max_frames",...
            "Mapping",{"lin",10,500}));
    end    
    methods
        function plugin = FreqConv
            plugin.input_buffer = dsp.AsyncBuffer;
            setup(plugin.input_buffer,[1,1]);
            plugin.output_buffer = dsp.AsyncBuffer;
            setup(plugin.output_buffer,[1,1]);
            plugin.overlap_buffer = zeros(plugin.buffer_size, plugin.channels);
            plugin.IR_frames_left = [];
            plugin.IR_frames_right = [];
            plugin.delay_fft_buffer_left = zeros(1024, 1000);
            plugin.delay_fft_buffer_right = zeros(1024, 1000);
        end
        function out = process(p, in)
            write(p.input_buffer, in);
            while p.input_buffer.NumUnreadSamples >= p.buffer_size
                input = read(p.input_buffer, p.buffer_size);
                w = p.write_head;
                p.delay_fft_buffer_left(:, w) = fft(input(:,1), p.NFFT);
                p.delay_fft_buffer_right(:, w) = fft(input(:,2), p.NFFT);
                w = w + 1;
                if w > 1000
                    w = 1;
                end
                p.write_head = w;
                overlap = zeros(p.buffer_size,2);
                r = w;

                % Dry signal%
                rev = input * (1-p.mix/100);

                % Wet signal
                rev = rev + p.overlap_buffer(1:p.buffer_size,:) * p.mix/100;
                n_frames = min(p.n_IR_frames, p.max_frames);
                for frame = 1:n_frames
                    in_left = p.delay_fft_buffer_left(:,r);
                    in_right = p.delay_fft_buffer_right(:,r);
                    IR_frame_real = real(p.IR_frames_left(1+(frame-1)*p.NFFT:frame*p.NFFT,1));
                    IR_frame_imag = imag(p.IR_frames_left(1+(frame-1)*p.NFFT:frame*p.NFFT,1));
                    res_left = FreqConvolute(in_left, IR_frame_real, IR_frame_imag, p.NFFT);
                    res_right = FreqConvolute(in_right, IR_frame_real, IR_frame_imag, p.NFFT);
                    rev(:,1) = rev(:,1) + res_left(1:p.buffer_size) / n_frames * p.mix/100;
                    rev(:,2) = rev(:,2) + res_right(1:p.buffer_size) / n_frames * p.mix/100;
                    overlap(:,1) = overlap(:,1) + res_left(p.buffer_size+1:p.buffer_size*2) / n_frames;
                    overlap(:,2) = overlap(:,2) + res_right(p.buffer_size+1:p.buffer_size*2) / n_frames;
                    if r == 1
                        r = 1000;
                    else
                        r = r - 1;
                    end
                end
                write(p.output_buffer,rev(:,1:2));
                p.overlap_buffer(1:p.buffer_size,:) = overlap(1:end,:);  
            end
            if p.output_buffer.NumUnreadSamples >= size(in,1)
                out = read(p.output_buffer,size(in,1));
            else
                out = zeros(size(in,1),2);
            end
        end
        
        function reset(p)
            IR_mat = coder.load("IR.mat"); 
            h = IR_mat.h_new;
            IR = RemoveTailBelowThreshold(h, -60);
            p.write_head = 1;
            [IR_frames_real, IR_frames_imag, ~] = GetUnisonPartitionedIRFrames(IR, p.NFFT, p.buffer_size);
            p.IR_frames_left = complex(IR_frames_real, IR_frames_imag);
            [IR_frames_real, IR_frames_imag, p.n_IR_frames] = GetUnisonPartitionedIRFrames(IR, p.NFFT, p.buffer_size);
            p.IR_frames_right = complex(IR_frames_real, IR_frames_imag);
            p.delay_fft_buffer_left = zeros(1024, 1000);
            p.delay_fft_buffer_right = zeros(1024, 1000);
            p.overlap_buffer = zeros(p.buffer_size, p.channels);
        end
    end
end