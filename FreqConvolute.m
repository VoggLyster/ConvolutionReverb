function [out] = FreqConvolute(x, ir_frame)
    len = numel(x);
    x = [x; zeros(numel(ir_frame)-len,1)];
    Y = ir_frame;
    NFFT = numel(Y);
    
    X = fft(x, NFFT);
    X_conv = X.*Y;
    out = real(ifft(X_conv,NFFT));
end
