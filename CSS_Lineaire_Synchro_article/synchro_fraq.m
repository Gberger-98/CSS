function [phi,lambda] = synchro_fraq(y,SF,B,alpha)
    
    Ts = 2^SF/B;
    M = 2*2^SF;
    te = 1/(2*B);

    z = [];
    for k=0:5
        tmp = conj(y((k+1)*M+1:(k+2)*M));
        z = [z, angle(sum((y(k*M+1:(k+1)*M).*tmp)))/(2*pi) ];
    end
    
    signe = sum(sign(z))/7;
    
    if abs(signe) ~= 1
    
        if signe <0
            phi = -mean(abs(z));
        else
            phi = mean(abs(z));
        end       
    else
        phi = mean(z);
    end

%     if phi >0.5
%         phi = -(1-phi);
%     end

%     if signe <0
%         phi = mean(abs(z));
%     else
%         phi = -mean(abs(z));
%     end
%     
    if abs(phi) >0.5
        phi = (1-phi);
    end


    y= y.*exp(-1j*2*pi*(-phi)*te/Ts*[0:length(y)-1]);
    
    tmp = [];
    
    for k=1:6
        a = fft(y(k*M+1:(k+1)*M),M);
        [~,faux_max] = max(abs(a));
        b = recherche_dichotomique(2*pi*(faux_max-2)/M,2*pi*(faux_max)/M,1e-4,y(k*M+1:(k+1)*M));
        tmp = [tmp,mod(b/(2*pi)*M,M)];
    end
    lambda =  alpha - floor((mean(tmp)-floor(mean(tmp)))*alpha);
    
end

