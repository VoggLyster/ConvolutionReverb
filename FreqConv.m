classdef FreqConv < audioPlugin
    properties
        IR_resample_percentage = 100;
        mix = 50;
        max_frames = 150;
    end
    properties (Access = private)
        IR_frames = [];
        n_IR_frames = 0;
        delay_buffer_size = 1000000;
        delay_buffer = [];
        channels = 2;
        buffer_size = 512;
        overlap_buffer = zeros(512,1);
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
        function out = process(p, in)
           frame_size = size(in,1);
           w = p.write_head;
           for n = 1:size(in,1)
              p.delay_buffer(w,:) = in(n,:);
              w = w + 1;
              if w > p.delay_buffer_size
                  w = 1;
              end
           end
           p.write_head = w;
           overlap = zeros(size(in));
           if w < frame_size
               r = w - frame_size + p.delay_buffer_size;
           else
               r = w - frame_size;
           end

           % Dry signal%
           out = p.delay_buffer(r:r+frame_size-1,:) * (1-p.mix/100);

           % Wet signal
           out = out + p.overlap_buffer(1:frame_size,:) * p.mix/100;
           n_frames = min(p.n_IR_frames, p.max_frames);
           for frame = 1:n_frames
              for channel = 1:p.channels
                  x = p.delay_buffer(r:r+frame_size-1,channel);
                  res = FreqConvolute(x, p.IR_frames(1:end,frame));
                  out(1:end,channel) = out(1:end,channel) + res(1:frame_size) / n_frames * p.mix/100;
                  overlap(1:end,channel) = overlap(1:end,channel) + res(frame_size+1:frame_size*2) / n_frames;
              end
              if r == 1
                r = p.delay_buffer_size - frame_size;
              else
                r = r - frame_size;
              end
           end
           p.overlap_buffer(1:frame_size,:) = overlap(1:end,:);  
        end
        
        function reset(p)
            IR_mat = coder.load("IR.mat"); 
            h = IR_mat.h_new;
            threshold = -60;
            threshold_mag = db2mag(threshold);
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
            p.IR_frames = GetIRFrames(IR, p.NFFT, p.buffer_size);
            p.n_IR_frames = size(p.IR_frames,2);
            p.delay_buffer = zeros(p.delay_buffer_size, p.channels);
            p.overlap_buffer = zeros(p.buffer_size, p.channels);
        end
        
        function set.IR_resample_percentage(p, val)
            p.IR_resample_percentage = val;
            p.reset();
        end
    end
end