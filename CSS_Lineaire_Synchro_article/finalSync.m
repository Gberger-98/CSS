function [deltaf_est,deltat_est,lsto] = finalSync(y,SF)
    M = 2^SF;
    [~ ,up] = max(abs(fft(y(4*M+1:5*M),M)));
    [~ ,down] = max(abs(fft(y(10*M+1:11*M),M)));
    
    
    deltaf_est = mod(round((up+down)),M)/2-1;
    
    if deltaf_est < M/4
        deltaf_est= mod(floor(deltaf_est),M);
        deltat_est = mod(floor(deltaf_est-up+1),M);
    else
        deltaf_est= mod(deltaf_est-M/2,M);
        deltat_est = mod(floor(deltaf_est-up+1),M);
    end
    
    Y= [];
    for m=[4,8]
        Y = [Y; fft(y(((m-1)*2^SF+1):(m*2^SF)),2^SF)];
    end
    M = 2^SF-deltat_est;
    Yavg = sum(Y); 
    
    tmp = exp(1j*2*pi*M/2^SF);
    if up ~=2^SF && up ~= 1 
        lsto = -real((conj(tmp)*Yavg(up+1)-tmp*Yavg(up-1))/(2*Yavg(up)-conj(tmp)*Yavg(up+1)-tmp*Yavg(up-1)));
    else
        if up == 1
            lsto = -real((conj(tmp)*Yavg(up+1)-tmp*Yavg(end))/(2*Yavg(up)-conj(tmp)*Yavg(up+1)-tmp*Yavg(end)));
        else
            lsto = -real((conj(tmp)*Yavg(1)-tmp*Yavg(up-1))/(2*Yavg(up)-conj(tmp)*Yavg(1)-tmp*Yavg(up-1)));
        end
    end
%     if abs(lsto)>0.5
%         lsto = -sign(lsto)*mod(lsto,0.5);
%     end
end

