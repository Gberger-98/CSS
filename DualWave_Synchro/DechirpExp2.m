function [Y_dechirp] = DechirpExp2(S,SF,B,alpha,tau,A,frac,pre)
    Ts = 1/(frac*B);
    T=frac*2^SF/(frac*B);
    t=[0:Ts:T-Ts]; 

    x_0 = exp(1j*2*pi*(A*(t+exp(-t/tau)*tau)-1*B*t/2));
    x_0_conj = conj(x_0);
    
    if alpha >1
        S = S(1:alpha/frac:end);
    else
        S = S(1:end);
    end
    N = floor(length(S)/(frac*2^SF));
    Y_dechirp = zeros(size(S));
    
    
    for i=1:N
        if (i==11 || i==12) && pre
            Y_dechirp(((i-1)*frac*2^SF+1):(i*frac*2^SF)) = S(((i-1)*frac*2^SF+1):(i*frac*2^SF)).*x_0;
        else
            Y_dechirp(((i-1)*frac*2^SF+1):(i*frac*2^SF)) = S(((i-1)*frac*2^SF+1):(i*frac*2^SF)).*x_0_conj;
        end
    end
    Y_dechirp((N*frac*2^SF):end)=S((N*frac*2^SF):end).*x_0_conj(1:length(S((N*frac*2^SF):end)));
    
end
