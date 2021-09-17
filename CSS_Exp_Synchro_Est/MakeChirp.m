function [Sy] = MakeChirp(SF,S,B,alpha)
    Ts = 1/(B*alpha);
    T=2^SF/B;
    t=[0:Ts:T-Ts]; 
    
    Sy = zeros(1,length(S)*length(t));
    
    x = exp(1j*2*pi*(B/(2*T)*t.^2-B*3/2*t));
    for i=1:length(S)
        nfold = alpha*(2^SF-S(i));
        Sy(alpha*((i-1)*2^SF)+1:alpha*(i*2^SF)) = circshift(x,nfold);
    end

end

