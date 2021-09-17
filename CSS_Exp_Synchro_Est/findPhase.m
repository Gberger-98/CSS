function [e] = findPhase(yc,SF,B,alpha,Nbpream)
    M=2^SF;
    T=M/B;
    Ts = 1/(alpha*B);
    t=[0:Ts:T-Ts];
    e = zeros(1,2^SF*alpha);
    for i=1:Nbpream
        sig_exp = yc(2^SF*alpha*(i-1)+1:alpha*2^SF*i);
        e = e + sig_exp;
    end
    e = e/Nbpream;
    
    
    
end

