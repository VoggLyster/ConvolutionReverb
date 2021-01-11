function [out, out_size] = FreqConvolute(input, ir_frame_real, ir_frame_imag, fft_frame_size)
    ir_frame = complex(ir_frame_real,ir_frame_imag);
    input_size = numel(input);
    input = [input; zeros(fft_frame_size-input_size,1)];
    Y = ir_frame;
    NFFT = fft_frame_size;
    
    X = fft(input, NFFT);
    X_conv = X.*Y;
    out = real(ifft(X_conv,NFFT));
    out_size = size(out,1);
end
