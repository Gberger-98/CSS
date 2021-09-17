function [Sy] = MakeChirp(SF,S,B,alpha)
    Ts = 2^SF/B; %Temps symbole
    te = 1/(alpha*B);
    t = (0:alpha*(2^SF)-1)*te; %Temps d'échantillonage 
    t2 = t.^2;
    
    Sy = zeros(1,length(S)*length(t));
    const1=B/(2*Ts);
    
    x = zeros(1,length(t));
    for i=1:length(S)
        nfold = alpha*(2^SF-S(i));
        x(1:nfold-1) = exp(1j*2*pi*(const1*t2(1:nfold-1)+B*(S(i)/2^SF-1/2)*t(1:nfold-1)));
        x(nfold:end) = exp(1j*2*pi*(const1*t2(nfold:end)+B*(S(i)/2^SF-3/2)*t(nfold:end)));
        Sy(alpha*((i-1)*2^SF)+1:alpha*(i*2^SF)) = x;
    end

end

