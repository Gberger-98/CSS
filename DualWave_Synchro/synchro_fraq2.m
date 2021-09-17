function [phi] = synchro_fraq2(y,SF,B)
    
    Ts = 2^SF/B;
    M = 2*2^SF;
    te = 1/(2*B);

    z = [];
    for k=0:2
        tmp = conj(y((k+1)*M+1:(k+2)*M));
        z = [z, angle(sum((y(k*M+1:(k+1)*M).*tmp)))/(2*pi) ];
    end
    
    signe = sum(sign(z))/3;
    
    if abs(signe) ~= 1
    
        if signe <0
            phi = -mean(abs(z));
        else
            phi = mean(abs(z));
        end       
    else
        phi = mean(z);
    end
    
    
end