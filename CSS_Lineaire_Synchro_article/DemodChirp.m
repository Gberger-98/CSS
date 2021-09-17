function [out] = DemodChirp(y,SF)
    
    
    N = floor(length(y)/2^SF);
    out=zeros(N,1);
    
    for i=1:N
        tmp = fft(y(((i-1)*2^SF+1):(i*2^SF)),2^SF);
        [~,k] = max(abs(tmp));
        out(i) = k-1;
    end

end

