function [deltaf_est,deltat_est] = synchro_int(y,SF,B)

    Ts = 2^SF/B;
    M = 2^SF;
    te = 1/B;
    
    
    %% Partie entière CFO STO
    [~ ,up] = max(abs(fft(y(4*M+1:5*M),M)));
    [~ ,down] = max(abs(fft(y(11*M+1:12*M),M)));
    
    
    deltaf_est = mod(round((up+down)),M)/2-1;
    
    if deltaf_est < M/4
        deltaf_est= mod(floor(deltaf_est),M);
        deltat_est = mod(floor(deltaf_est-up+1),M);
    else
        deltaf_est= mod(deltaf_est-M/2,M);
        deltat_est = mod(floor(deltaf_est-up+1),M);
    end

end

