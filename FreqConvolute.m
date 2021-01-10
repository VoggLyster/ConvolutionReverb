function [out, out_size] = FreqConvolute(input, input_size, ir_frame_real, ir_frame_imag, frame_size)
    ir_frame = complex(ir_frame_real,ir_frame_imag);
    input = [input; zeros(frame_size-input_size,1)];
    Y = ir_frame;
    NFFT = frame_size;
    
    X = fft(input, NFFT);
    X_conv = X.*Y;
    out = real(ifft(X_conv,NFFT));
    out_size = size(out,1);
end
