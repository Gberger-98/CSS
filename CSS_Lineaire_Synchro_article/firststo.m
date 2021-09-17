function [lsto] = firststo(y,SF)
    Y= [];
    for m=1:3
        Y = [Y; fft(y(((m-1)*2^SF+1):(m*2^SF)),2^SF)];
    end
    
    [~,Sup] = max(abs(Y(1,:)));
    M = 2^SF-Sup+1;
    
    Yavg = sum(Y); 
    
    tmp = exp(1j*2*pi*M/2^SF);
    if Sup ~=2^SF && Sup ~= 1 
        lsto = -real((conj(tmp)*Yavg(Sup+1)-tmp*Yavg(Sup-1))/(2*Yavg(Sup)-conj(tmp)*Yavg(Sup+1)-tmp*Yavg(Sup-1)));
    else
        if Sup == 1
            lsto = -real((conj(tmp)*Yavg(Sup+1)-tmp*Yavg(end))/(2*Yavg(Sup)-conj(tmp)*Yavg(Sup+1)-tmp*Yavg(end)));
        else
            lsto = -real((conj(tmp)*Yavg(1)-tmp*Yavg(Sup-1))/(2*Yavg(Sup)-conj(tmp)*Yavg(1)-tmp*Yavg(Sup-1)));
        end
    end
%     if abs(lsto)>0.5
%         lsto = -sign(lsto)*mod(lsto,0.5);
%     end
    
end

