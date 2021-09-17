function [Sy] = MakeChirpExp2(SF,S,B,alpha,tau,A)
    Ts = 1/(B*alpha);
    T=2^SF/B;
    t=[0:Ts:T-Ts]; 
    
    Sy = zeros(1,length(S)*length(t));
    
%     x = exp(1j*2*pi*(A*(t+exp(-t/tau)*tau)+B*(0/2^SF-3/2)*t));
    x = zeros(1,length(t));
    for i=1:length(S)
        nfold = alpha*(2^SF-S(i));
        x(1:nfold) = exp(1j*2*pi*(A*(t(1:nfold)+exp(-t(1:nfold)/tau)*tau)+B*(S(i)/2^SF-1/2)*t(1:nfold)));
        x(nfold:end) = exp(1j*2*pi*(A*(t(nfold:end)+exp(-t(nfold:end)/tau)*tau)+B*(S(i)/2^SF-3/2)*t(nfold:end)));
%         if S(i) == 0
%             x = exp(1j*2*pi*(A*(t+exp(-t/tau)*tau)+B*(0/2^SF-3/2)*t));
%         end
        Sy(alpha*((i-1)*2^SF)+1:alpha*(i*2^SF)) = x;%circshift(x,nfold);
    end

end

