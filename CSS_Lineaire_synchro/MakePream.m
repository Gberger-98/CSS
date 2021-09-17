function [pream] = MakePream(SF,B,alpha)
    seq = [zeros(8,1);4;4];
    pream=MakeChirp(SF,seq,B,alpha);
    downchirp = conj(pream(1:alpha*2^SF));
    pream = [pream, downchirp,downchirp];

end

