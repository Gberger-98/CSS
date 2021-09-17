function [Y_dechirp] = Dechirp(S,SF,B,alpha,pre)
    M=2^SF;
    T=M/B;
    Ts = 1/B;
    t=[0:Ts:T-Ts];
    x_0 = exp(1j*2*pi*(B/(2*T)*t.^2-B*3/2*t));
%     x_0 = exp(1j*2*pi*(B/(2*T)*t.^2));
    x_0_conj = conj(x_0);
    
    if alpha >1
        S = S(1:alpha:end);
    else
        S = S(1:end);
    end
    N = floor(length(S)/2^SF);
    Y_dechirp = zeros(size(S));
    
    
    for i=1:N
        if (i==11 || i==12) && pre
            Y_dechirp(((i-1)*2^SF+1):(i*2^SF)) = S(((i-1)*2^SF+1):(i*2^SF)).*x_0;
        else
            Y_dechirp(((i-1)*2^SF+1):(i*2^SF)) = S(((i-1)*2^SF+1):(i*2^SF)).*x_0_conj;
        end
    end
    Y_dechirp((N*2^SF):end)=S((N*2^SF):end).*x_0_conj(1:length(S((N*2^SF):end)));
    
end
