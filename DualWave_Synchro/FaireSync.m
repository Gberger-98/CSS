function [start] = FaireSync(yc,pream,alpha,SF)
    
    LengthSeqGold = length(pream);
    M = 2^SF;
    n_corr=alpha*M*10+LengthSeqGold;
    energyRefSig = pream*pream';
    refSig = pream;
    absCplxBuffer = abs(yc(1:2*n_corr));
    
    NumerateurSynch = filter(fliplr(pream),1,absCplxBuffer);
    Etestsig = filter(ones(1,length(pream)),1,absCplxBuffer.^2);
    DenominateurSynch = sqrt(energyRefSig*Etestsig);
    corrVal = NumerateurSynch./DenominateurSynch;
    [valMaxSynch,indMaxSynch] = max(abs(corrVal));
    start = indMaxSynch-length(pream);
    
    
end

