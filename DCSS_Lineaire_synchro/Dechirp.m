function [Y_dechirp] = Dechirp(S,SF,B,alpha,frac)
    Ts = 2^SF/B;
    t = (0:frac*2^SF-1)*1/(frac*B);
    
    x_0 = exp(1j*2*pi*(B/(2*Ts)*t.^2-B/2*t));
    x_0_conj = conj(x_0);
    
    if alpha >1
        S = S(1:alpha/frac:end);
    else
        S = S(1:end);
    end
    N = floor(length(S)/(frac*2^SF));
    Y_dechirp = zeros(size(S));
    
    
    for i=1:N
        Y_dechirp(((i-1)*frac*2^SF+1):(i*frac*2^SF)) = S(((i-1)*frac*2^SF+1):(i*frac*2^SF)).*x_0_conj;
    end
    Y_dechirp((N*frac*2^SF):end)=S((N*frac*2^SF):end).*x_0_conj(1:length(S((N*frac*2^SF):end)));
    
end
