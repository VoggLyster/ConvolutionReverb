function [ir_frames_left, ir_frames_right, n_frames] = GetUnisonPartitionedIRFrames(IR, NFFT, buffer_size)
    IR_stereo = zeros(size(IR,1),2);
    if size(IR,2) ~= 2
       IR_stereo(:,1) = IR(:,1); 
       IR_stereo(:,2) = IR(:,1);
    else 
       IR_stereo = IR(:,1:2);
    end
    n_frames = ceil(size(IR_stereo,1)/buffer_size);
    IR_frames_left = zeros(buffer_size, n_frames);
    IR_frames_right = zeros(buffer_size, n_frames);
    IR_fft_frames_left = complex(zeros(NFFT*n_frames, 1));
    IR_fft_frames_right = complex(zeros(NFFT*n_frames, 1));
    IR_length_ceil = n_frames * buffer_size;
    padding_needed = IR_length_ceil - size(IR_stereo,1);
    IR_padded = vertcat(IR_stereo, zeros(padding_needed,2));
    idx_start = 1;
    idx_end = buffer_size;
    for frame = 1:n_frames
        IR_frames_left(1:end,frame) = IR_padded(idx_start:idx_end,1);
        IR_frames_right(1:end,frame) = IR_padded(idx_start:idx_end,2);
        idx_start = idx_end+1;
        idx_end = idx_end + buffer_size;
    end
    for frame = 1:n_frames
        IR_frame_len_left = length(IR_frames_left(1:end,frame));
        IR_frame_len_right = length(IR_frames_right(1:end,frame));
        IR_frame_padded_left = [IR_frames_left(1:end,frame); zeros(NFFT - IR_frame_len_left, 1)];
        IR_frame_padded_right = [IR_frames_right(1:end,frame); zeros(NFFT - IR_frame_len_right, 1)];
        IR_frame_left = fft(IR_frame_padded_left,NFFT);
        IR_frame_right = fft(IR_frame_padded_right,NFFT);
        IR_fft_frames_left(1+(frame-1)*NFFT:frame*NFFT,1) = IR_frame_left;
        IR_fft_frames_right(1+(frame-1)*NFFT:frame*NFFT,1) = IR_frame_right;
    end
    ir_frames_left = IR_fft_frames_left;
    ir_frames_right = IR_fft_frames_right;
end