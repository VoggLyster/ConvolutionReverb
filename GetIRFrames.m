function [ir_frames] = GetIRFrames(IR, NFFT, buffer_size)
    n_frames = ceil(size(IR,1)/buffer_size);
    IR_frames_left = zeros(buffer_size, n_frames);
    IR_fft_frames_left = zeros(NFFT, n_frames);
    IR_length_ceil = n_frames * buffer_size;
    padding_needed = IR_length_ceil - size(IR,1);
    IR_padded = vertcat(IR, zeros(padding_needed,1));
    idx_start = 1;
    idx_end = buffer_size;
    for frame = 1:n_frames
        IR_frames_left(1:end,frame) = IR_padded(idx_start:idx_end,1);
        idx_start = idx_end+1;
        idx_end = idx_end + buffer_size;
    end
    for frame = 1:n_frames
        IR_frame_len = length(IR_frames_left(1:end,frame));
        IR_frame_padded = [IR_frames_left(1:end,frame); zeros(NFFT - IR_frame_len, 1)];
        IR_frame = fft(IR_frame_padded,NFFT);
        IR_fft_frames_left(1:end,frame) = IR_frame;
    end
    ir_frames = IR_fft_frames_left;
end