classdef FreqConv < audioPlugin
    properties
        IR_resample_percentage = 100;
        mix = 50;
        max_frames = 150;
    end
    properties (Access = private)
        input_buffer;
        output_buffer;
        overlap_buffer;
        delay_buffer;
        IR_frames_left;
        IR_frames_right;
        n_IR_frames = 0;
        delay_buffer_size = 1000000;
        channels = 2;
        buffer_size = 512;
        write_head = 1;
        NFFT = 0;
    end
    properties (Constant)
        PluginInterface = audioPluginInterface( ...
            audioPluginParameter("IR_resample_percentage", ...
            "Mapping",{"lin",20,200}),...
            audioPluginParameter("mix",...
            "Mapping",{"lin",0,100}),...
            audioPluginParameter("max_frames",...
            "Mapping",{"lin",10,1000}));
    end    
    methods
        function plugin = FreqConv
           plugin.input_buffer = dsp.AsyncBuffer;
           setup(plugin.input_buffer,[1,1]);
           plugin.output_buffer = dsp.AsyncBuffer;
           setup(plugin.output_buffer,[1,1]);
           plugin.delay_buffer = zeros(plugin.delay_buffer_size, plugin.channels);
           plugin.overlap_buffer = zeros(plugin.buffer_size, plugin.channels);
        end
        function out = process(p, in)
           write(p.input_buffer, in);
           while p.input_buffer.NumUnreadSamples >= p.buffer_size
               w = p.write_head;
               p.delay_buffer(w:w+p.buffer_size-1,:) = read(p.input_buffer, p.buffer_size);
               w = w + p.buffer_size;
               if w > p.delay_buffer_size
                  w = 1;
               end
               p.write_head = w;
               overlap = zeros(p.buffer_size,2);
               if w < p.buffer_size
                   r = w - p.buffer_size + p.delay_buffer_size;
               else
                   r = w - p.buffer_size;
               end

               % Dry signal%
               rev = p.delay_buffer(r:r+p.buffer_size-1,:) * (1-p.mix/100);

               % Wet signal
               rev = rev + p.overlap_buffer(1:p.buffer_size,:) * p.mix/100;
               n_frames = min(p.n_IR_frames, p.max_frames);
               for frame = 1:n_frames
                      x = p.delay_buffer(r:r+p.buffer_size-1,:);
                      res_left = FreqConvolute(x(:,1), p.IR_frames_left(1+(frame-1)*p.NFFT:frame*p.NFFT,1));
                      res_right = FreqConvolute(x(:,2), p.IR_frames_right(1+(frame-1)*p.NFFT:frame*p.NFFT,1));
                      rev(:,1) = rev(:,1) + res_left(1:p.buffer_size) / n_frames * p.mix/100;
                      rev(:,2) = rev(:,2) + res_right(1:p.buffer_size) / n_frames * p.mix/100;
                      overlap(:,1) = overlap(:,1) + res_left(p.buffer_size+1:p.buffer_size*2) / n_frames;
                      overlap(:,2) = overlap(:,2) + res_right(p.buffer_size+1:p.buffer_size*2) / n_frames;
                  if r == 1
                    r = p.delay_buffer_size - p.buffer_size;
                  else
                    r = r - p.buffer_size;
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
            threshold = -60;
            threshold_mag = db2mag(threshold);
            IR = [];
            for i = numel(h):-1:1
               if abs(h(i)) > threshold_mag
                  IR = h(1:i);
                  break 
               end
            end
            p.NFFT = 2^nextpow2(p.buffer_size + 1);
            p.write_head = 1;
            if(p.IR_resample_percentage ~= 1)
                IR = resample(IR,round(p.IR_resample_percentage),100);
            end
            [p.IR_frames_left, p.IR_frames_right, p.n_IR_frames] = GetUnisonPartitionedIRFrames(IR, p.NFFT, p.buffer_size);
            p.delay_buffer = zeros(p.delay_buffer_size, p.channels);
            p.overlap_buffer = zeros(p.buffer_size, p.channels);
        end
        
        function set.IR_resample_percentage(p, val)
            p.IR_resample_percentage = val;
            p.reset();
        end
    end
end