function [Y_dechirp] = DechirpExp3(S,SF,B,alpha,a)
    Ts = 1/B;
    T=2^SF/B;
    t=[0:Ts:T-Ts]; 
    a = a(1:alpha:end);
    x_0 = exp(1j*a);
    x_0_conj = conj(x_0);
    
    if alpha >1
        S = S(1:alpha:end);
    else
        S = S(1:end);
    end
    N = floor(length(S)/2^SF);
    Y_dechirp = zeros(size(S));
    
    
    for i=1:N
        Y_dechirp(((i-1)*2^SF+1):(i*2^SF)) = S(((i-1)*2^SF+1):(i*2^SF)).*x_0_conj;
    end
    Y_dechirp((N*2^SF):end)=S((N*2^SF):end).*x_0_conj(1:length(S((N*2^SF):end)));
    
end

