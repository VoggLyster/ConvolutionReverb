function [ir_frames_real, ir_frames_imag, n_frames, output_size] = GetUnisonPartitionedIRFrames(IR, IR_size, frame_size, buffer_size, ir_frames_real, ir_frames_imag)
    n_frames = ceil(IR_size/buffer_size);
    IR_frames = zeros(buffer_size, n_frames);
    IR_fft_frames = complex(zeros(frame_size*n_frames, 1));
    IR_length_ceil = n_frames * buffer_size;
    padding_needed = IR_length_ceil - IR_size;
    IR_padded = vertcat(IR(1:IR_size), zeros(padding_needed,1));
    idx_start = 1;
    idx_end = buffer_size;
    for frame = 1:n_frames
        IR_frames(1:end,frame) = IR_padded(idx_start:idx_end,1);
        idx_start = idx_end+1;
        idx_end = idx_end + buffer_size;
    end
    for frame = 1:n_frames
        this_frame_size = length(IR_frames(1:end,frame));
        IR_frame_padded = [IR_frames(1:end,frame); zeros(frame_size - this_frame_size, 1)];
        IR_frame = fft(IR_frame_padded,frame_size);
        IR_fft_frames(1+(frame-1)*frame_size:frame*frame_size,1) = IR_frame;
    end
    output_size = n_frames * frame_size;
    ir_frames_real(1:output_size) = real(IR_fft_frames);
    ir_frames_imag(1:output_size) = imag(IR_fft_frames);
    
end