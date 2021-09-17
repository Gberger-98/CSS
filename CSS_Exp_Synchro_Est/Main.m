clc;
% close all;
clear;

%% Paramètre

B = 125e3; %Largeur de bande
% SF = [12,11,10,9,8,7]; %Nombre de bits par symbole
% SF = [9,8,7];
SF = 7;


alpha = 10;
B_E=[-1 0.4]; %Paramètre de non-linéarité dans ]0,1[ ou -1 pour chirp linéaire
% B_E = [-1, 0.2 0.5 0.7];
% B_E = 0.6;

Nbbit = 1e6;

N = 100; %Nombre de symbole



Nbpream = 10;

eb=1;
% eb_n0_dB = 100; %Rapport signal sur bruit
eb_n0_dB = -10:0;
eb_n0    = 10 .^ (eb_n0_dB/10);         % Liste des Eb/N0
sigma2   = 1 ./eb_n0; % Variance du bruit complex


TEB = zeros(length(eb_n0),length(SF),length(B_E));   % Tableau des TEB (resultats)
PER = zeros(length(eb_n0),length(SF),length(B_E));
c = [];


for l=1:length(B_E)
    for k=1:length(SF)
        M=2^SF(k);
        T=M/B;
        g = hann(round(SF(k)));
        g = g/sum(g);
        if B_E(l) ~= -1
            tau = -T/log(1-B_E(l));
            A = B/B_E(l);
        else
            tau = 0;
        end
        c = [c, B_E(l)];
        for i = 1:length(eb_n0)
            error_cnt   = 0;
            bit_cnt     = 0;
            pck_cnt     = 0;
            err_pck_cnt = 0;
            if B_E(l) == -1
                pream = MakeChirp(SF(k),zeros(Nbpream,1),B,alpha);
            else
                pream = MakeChirpExp2(SF(k),zeros(Nbpream,1),B,alpha,tau,A);
            end
            while bit_cnt < Nbbit
                %% Generation
                bitsM = randi([0,1],[N,SF(k)]);
                numsM = bi2de(bitsM);

                
                %% Modulation
                if B_E(l) == -1
                    S = MakeChirp(SF(k),numsM,B,alpha);
                else
                    S = MakeChirpExp2(SF(k),numsM,B,alpha,tau,A);
                end
                
                %% Canal
                S = [pream,S];
                
                
                wl = sqrt(sigma2(i)/2)*(randn(size(S))+ 1j*randn(size(S)));
                yc = S + wl;
                
                
                
                %% Demodulation
                if B_E(l) ~= -1
                    [etmp] = findPhase(yc,SF(k),B,alpha,Nbpream); % Estimation du downchirp
                    if exist('e','var')== 1
                        e = (e+ etmp)/2;
                    else
                        e = etmp;
                    end
                    ph = phase(e);
                    a = conv(ph, g, 'same'); % Moyennage
                    
                    y = DechirpExp3(yc,SF(k),B,alpha,a);
                else
                    y = Dechirp(yc,SF(k),B,alpha,0);
                end
                nums_est = DemodChirp(y,SF(k));
                
                
                Offset = round(mean(nums_est(1:Nbpream)));
                SymbEstCorr = mod((nums_est-Offset),2^SF(k));
                bitsM_est = de2bi(SymbEstCorr(Nbpream+1:end),SF(k));
                %             tmp = sum(nums_est ~= numsM,"all");
                tmp = sum(bitsM_est(1:end,:) ~= bitsM(1:end,:),"all");
                
                
                bit_cnt = bit_cnt + N*SF(k);
                
                error_cnt =  error_cnt + tmp;
                
                pck_cnt     = pck_cnt + 1;
                if tmp >0
                    err_pck_cnt = err_pck_cnt + 1;
                end
                
                
            end
            TEB(i,k,l) = error_cnt/bit_cnt;
            PER(i,k,l) = err_pck_cnt/pck_cnt;
            
            fprintf("TEB = %1.2e à SNR = %2.1fdB et SF = %i et B_E= %1.2f\n",TEB(i,k,l),eb_n0_dB(i),SF(k), (B_E(l)~=-1)*B_E(l) )
            clear e;
            if TEB(i,k,l) == 0 %&& B_E(l) == -1
                break
            end
            
        end
    end
end

%% Affichage
figure
semilogy(eb_n0_dB,TEB(:,1,1))
hold
for m=2:length(B_E)
    semilogy(eb_n0_dB,TEB(:,1,m),'--o')
end
styl = ["--","--s","--*","--x","--d","--p","--^","-->","--<","--h"];
for i=2:length(SF)
    for m=1:length(B_E)
        if m ==1
            semilogy(eb_n0_dB,TEB(:,i,m))
        else
            semilogy(eb_n0_dB,TEB(:,i,m),styl(i))
        end
    end
end

%% Legend
Legend = strings(length(SF)*length(B_E),1);

for i=1:length(SF)
    for m=1:length(B_E)
        if m==1
            Legend((i-1)*length(B_E)+m) = sprintf("SF = %i ",SF(i)); 
        else
            Legend((i-1)*length(B_E)+m) = sprintf("B_E = %1.2f",B_E(m));
        end
    end
end
legend(Legend);
grid ON
ylabel("BER")
xlabel("SNR (dB)")

%% PER
figure
semilogy(eb_n0_dB,PER(:,1,1))
hold
for m=2:length(B_E)
    semilogy(eb_n0_dB,PER(:,1,m),'--o')
end
styl = ["--","--s","--*","--x","--d","--p","--^","-->","--<","--h"];
for i=2:length(SF)
    for m=1:length(B_E)
        if m ==1
            semilogy(eb_n0_dB,PER(:,i,m))
        else
            semilogy(eb_n0_dB,PER(:,i,m),styl(i))
        end
    end
end

%% Legend
Legend = strings(length(SF)*length(B_E),1);

for i=1:length(SF)
    for m=1:length(B_E)
        if m==1
            Legend((i-1)*length(B_E)+m) = sprintf("SF = %i ",SF(i)); 
        else
            Legend((i-1)*length(B_E)+m) = sprintf("B_E = %1.2f",B_E(m));
        end
    end
end
legend(Legend);
grid ON
ylabel("PER")
xlabel("SNR (dB)")


