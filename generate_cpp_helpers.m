max_IR_size = 2^20;
buffer_size = 2048;
IR_mat = coder.load("IR.mat"); 
h = IR_mat.h_new;
IR = RemoveTailBelowThreshold(h, -60);
fft_frame_size = 2^nextpow2(buffer_size + 1);
[IR_frames_real, IR_frames_imag, n_IR_frames] = GetUnisonPartitionedIRFrames(IR, fft_frame_size, buffer_size);
x = rand(buffer_size,1);
FreqConvolute(x, IR_frames_real(1:frame_size), IR_frames_imag(1:frame_size), fft_frame_size);