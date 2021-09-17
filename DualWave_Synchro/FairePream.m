function [pream] = FairePream(SF,B,alpha)
    
    M=2^SF;
    Ts = M/B;
    Te=Ts/(alpha*M);
    G = 5; % param�tre permettant de r�gler la localisation en fr�quence de l'�nergie de la DSP dans B
    Dbp = B/G;
    Fse = 1/(Dbp*Te);
    g = ones(1,Fse);
    Nbp = floor(10*M/G);
%     Nbp = 1024;
    name = sprintf('gold%i.mat',2^(SF+1));
%     name = 'gold1024.mat';
    goldSeq256 = load(name);
    sp(1,:) = filter(g,1,upsample(goldSeq256.X(:,1).',Fse));
%     sp(2,:) = filter(g,1,upsample(goldSeq256.X(:,2).',Fse));
%     sp(3,:) = filter(g,1,upsample(goldSeq256.X(:,3).',Fse));
%     sp(4,:) = filter(g,1,upsample(goldSeq256.X(:,4).',Fse));
%     sp(5,:) = filter(g,1,upsample(goldSeq256.X(:,5).',Fse));
    pream = [sp(1,:)];
    
    

end

