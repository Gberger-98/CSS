function [out] = DemodChirp(y,SF,flag)
    
    M = 2^SF;
    N = floor(length(y)/2^SF);
    out=zeros(N,1);
    epsilon = 10^(-4);
    
    for i=1:N
        sig = y(((i-1)*2^SF+1):(i*2^SF));
        tmp = fft(sig,2^SF);
        [~,k] = max(abs(tmp));
        if flag == 1
            a = 2*pi*(k-2)/M;
            b = 2*pi*k/M;
            [pos] = recherche_dichotomique(a, b, epsilon, sig);
            out(i) = pos/(2*pi)*M;
        else
            out(i)=k-1;
        end
    end
end

