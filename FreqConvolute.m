function [out, out_size] = FreqConvolute(input, ir_frame_real, ir_frame_imag, fft_frame_size)
    ir_frame = complex(ir_frame_real,ir_frame_imag);
    Y = ir_frame;
    NFFT = fft_frame_size;
    X_conv = input.*Y;
    out = real(ifft(X_conv,NFFT));
    out_size = size(out,1);
end
