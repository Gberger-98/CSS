clc;
% close all;
clear;

%% Paramètre


% SF = [12,11,10,9,8,7]; %Nombre de bits par symbole
% SF = [9,8,7];
SF = 7;



alpha = 10;
B = 125e3; %Largeur de bande
te = 1/(alpha*B);
N = 100; %Nombre de symbole
B_E=[-2 -1 0.3]; %Paramètre de non-linéarité dans ]0,1[ ou -1 pour chirp linéaire avec synchro, -2 pour chirp linéaire avec synchro parfaite
% B_E = -1;

eb=1;
% eb_n0_dB = 100; %Rapport signal su bruit
eb_n0_dB = -10:0;
eb_n0    = 10 .^ (eb_n0_dB/10);         % Liste des Eb/N0
sigma2   = 1 ./eb_n0; % Variance du bruit complex


TEB = zeros(length(eb_n0),length(SF),length(B_E));   % Tableau des TEB (resultats)
PER = zeros(length(eb_n0),length(SF),length(B_E));
c = [];
%% défaut

maxcfo = B/4.2;
nbmax = 1e6;

for l=1:length(B_E)
    for k=1:length(SF)
        M=2^SF(k);
        Ts = 2^SF(k)/B;
        T=M/B;
        if B_E(l) ~= -1 && B_E(l) ~= -2
            tau = -T/log(1-B_E(l));
            A = B/B_E(l);
        else
            tau = 0;
        end
        c = [c, B_E(l)];
        
        pream = FairePream(SF(k),B,alpha);
        if B_E(l) == -1 || B_E(l) == -2 
            rawchiprs = MakeChirp(SF(k),[0 0 0],B,alpha);
        else
            rawchiprs = MakeChirpExp2(SF(k),[0 0 0],B,alpha,tau,A);
        end
        
        for i = 1:length(eb_n0)
            error_cnt = 0;
            bit_cnt   = 0;
            pck_cnt   = 0;
            err_pck   = 0;
            while bit_cnt < nbmax && error_cnt < nbmax/8 
                %% Generation
                bitsM = randi([0,1],[N,SF(k)]);

                numsM = bi2de(bitsM);
                DnumsM = Diffp(numsM,SF(k));

                %% Modulation

                if B_E(l) == -1 || B_E(l) == -2
                    S = MakeChirp(SF(k),DnumsM,B,alpha);
                else
                    S = MakeChirpExp2(SF(k),DnumsM,B,alpha,tau,A);
                end
                S = [pream,rawchiprs,S];
                %% Canal


                deltat = randi([0,alpha*2^SF(k)/4-1],1);
                deltaf = randi(round(maxcfo));
                if B_E(l) ~= -2
                    S = [zeros(1,deltat),S];
                    wl = sqrt(sigma2(i)/2)*(randn(size(S))+ 1j*randn(size(S)));
                
                    yc = S.*exp(1j*2*pi*deltaf*te*[0:length(S)-1])+ wl;
                else
                    wl = sqrt(sigma2(i)/2)*(randn(size(S))+ 1j*randn(size(S)));
                    yc = S+ wl;
                end

                %% Demod
                
                if B_E(l) == -2
                    y3 = Dechirp(yc(length(pream)+1:end),SF(k),B,alpha,1);
                else
                    %% Synchro temps
                    ns = FaireSync(yc,pream,alpha,SF(k));
                    y2 = yc(ns+length(pream)+1:end);
                    if B_E(l) == -1
                        y3 = Dechirp(y2,SF(k),B,alpha,2);
                    else
                        y3 = DechirpExp2(y2,SF(k),B,alpha,tau,A,2,0);
                    end
                    phi = synchro_fraq2(y3,SF(k),B);
                    ytmp = y2.*exp(-1j*2*pi*(-phi)*te/Ts*[0:length(y2)-1]);
                    if B_E(l) == -1
                        y3 = Dechirp(ytmp,SF(k),B,alpha,1);
                    else
                        y3 = DechirpExp2(ytmp,SF(k),B,alpha,tau,A,1,0);
                    end
                end
                y3 = y3(3*2^SF(k)+1:end);

                if B_E(l) == -2
                    dnums_est = DemodChirp(y3,SF(k),0);
                else
                    dnums_est = DemodChirp(y3,SF(k),1);
                end
                
                
                nums_est2 = DeDiffp(dnums_est,SF(k));
                if length(nums_est2)>N
                    nums_est=nums_est2(1:N);
                else
                    nums_est=nums_est2;
                end
                bitsM_est = de2bi(mod(nums_est,2^SF(k)),SF(k));

                tmp = sum(bitsM_est(2:end,:) ~= bitsM(2:length(bitsM_est),:),"all");

                if tmp > 0
                    err_pck = err_pck + 1;
                end
                pck_cnt = pck_cnt+1;

                bit_cnt = bit_cnt + N*SF(k);
                error_cnt =  error_cnt + tmp;


            end
            TEB(i,k,l) = error_cnt/bit_cnt;
            PER(i,k,l) = err_pck/pck_cnt;
            
            fprintf("TEB = %1.2e à SNR = %2.1fdB et SF = %i et B_E= %1.2f\n",TEB(i,k,l),eb_n0_dB(i),SF(k),B_E(l))
            
            if TEB(i,k,l) == 0
                break
            end
        end
    end
end
%% BER
figure
semilogy(eb_n0_dB,TEB(:,1,1))
hold
semilogy(eb_n0_dB,TEB(:,1,2),'-o')
for m=3:length(B_E)
    semilogy(eb_n0_dB,TEB(:,1,m),'--d')
end
styl = ["--","--s","--*","--x","--d","--p","--^","-->","--<","--h"];
for i=2:length(SF)
    for m=1:length(B_E)
        if m ==1
            semilogy(eb_n0_dB,TEB(:,i,m))
        else
            if B_E(m) == -1
                semilogy(eb_n0_dB,TEB(:,i,m),"-o")
            else
                semilogy(eb_n0_dB,TEB(:,i,m),styl(i))
            end
        end
    end
end


Legend = strings(length(SF)*length(B_E),1);

for i=1:length(SF)
    for m=1:length(B_E)
        if m==1
            Legend((i-1)*length(B_E)+m) = sprintf("SF = %i with perfect synchronisation",SF(i)); 
        else
            if B_E(m) == -1
                Legend((i-1)*length(B_E)+m) = "linear";
            else
                Legend((i-1)*length(B_E)+m) = sprintf("B_E = %1.2f",B_E(m));
            end
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
semilogy(eb_n0_dB,PER(:,1,2),'-o')
for m=3:length(B_E)
    semilogy(eb_n0_dB,PER(:,1,m),'--d')
end
styl = ["--","--s","--*","--x","--d","--p","--^","-->","--<","--h"];
for i=2:length(SF)
    for m=1:length(B_E)
        if m ==1
            semilogy(eb_n0_dB,PER(:,i,m))
        else
            if B_E(l) == -1
                semilogy(eb_n0_dB,TEB(:,i,m),"-o")
            else
                semilogy(eb_n0_dB,TEB(:,i,m),styl(i))
            end
        end
    end
end


Legend = strings(length(SF)*length(B_E),1);

for i=1:length(SF)
    for m=1:length(B_E)
        if m==1
            Legend((i-1)*length(B_E)+m) = sprintf("SF = %i with perfect synchronisation ",SF(i)); 
        else
            if B_E(m) == -1
                Legend((i-1)*length(B_E)+m) = "linear";
            else
                Legend((i-1)*length(B_E)+m) = sprintf("B_E = %1.2f",B_E(m));
            end
        end
    end
end
legend(Legend);
grid ON
ylabel("PER")
xlabel("SNR (dB)")

