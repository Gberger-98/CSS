function [ns] = synchro_int(y,SF,B,alpha)

    T = 2^SF/B;
    M = 2^SF;
    Ts = 1/(alpha*B);
    t=[0:Ts:T-Ts];
    x_0 = exp(1j*2*pi*(B/(2*Ts)*t.^2-B/2*t));
    
    
    ns = binary_search(y,1,M,10,alpha,SF,x_0);

end

