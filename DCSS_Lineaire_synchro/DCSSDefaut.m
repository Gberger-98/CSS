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
N = 1000; %Nombre de symbole


eb=1;
eb_n0_dB = -5; %Rapport signal su bruit
% eb_n0_dB = 100;
eb_n0    = 10 .^ (eb_n0_dB/10);         % Liste des Eb/N0
sigma2   = 1 ./eb_n0; % Variance du bruit complex


TEB = zeros(length(eb_n0),length(SF));   % Tableau des TEB (resultats)
TEP = zeros(length(eb_n0),length(SF));
%% défaut

maxcfo = B/4.2;


for k=1:length(SF)
    pream1 = MakeChirp(SF(k),[zeros(1,10) 4 8],B,alpha);
    pream = [pream1 zeros(1,alpha*2^SF(k)/2)];
    Ts = 2^SF(k)/B;
    for i = 1:length(eb_n0)
        error_cnt = 0;
        bit_cnt   = 0;
        pck_cnt   = 0;
        err_pck   = 0;
        while bit_cnt < 1e5
            %% Generation
            bitsM = randi([0,1],[N,SF(k)]);

            numsM = bi2de(bitsM);
            DnumsM = Diffp(numsM,SF(k));

            %% Modulation
            
            S = MakeChirp(SF(k),DnumsM,B,alpha);
            S = [pream,S];
            %% Canal
            

            deltat = randi([0,alpha*(2^SF(k))/2-1],1);
            deltaf = randi(round(maxcfo));
            if deltat>=0
                S = [zeros(1,alpha*2^SF(k)+deltat),S];
            end
            wl = sqrt(sigma2(i)/2)*(randn(size(S))+ 1j*randn(size(S)));
            yc = S.*exp(1j*2*pi*deltaf*te*[0:length(S)-1])+ wl;

            %% Demodulation
           
            [ns] = synchro_int(yc,SF(k),B,alpha);
            y = yc(alpha*ns:end);
            y = Dechirp(y,SF(k),B,alpha,2);
            [phi,lambda] = synchro_fraq(y,SF(k),B,alpha);
            
            if ns >2^SF(k)
                deb = alpha*ns+lambda;
            else
                deb = alpha*2^SF(k)+lambda;
            end
            y3 = yc(deb:end).*exp(-1j*2*pi*(-phi)*te/Ts*[0:length(yc(deb:end))-1]);
            
            y4 = Dechirp(y3(1:length(pream)-alpha*2^SF(k)/2),SF(k),B,alpha,1);
            deltat2 = deb - alpha*2^SF(k);
            
%              fprintf("delta f =  %i, est = %i, dif = %i \n",deltaf*Ts,deltaf_est-phi,deltaf*Ts - (deltaf_est-phi))
%             fprintf("delta t = %i, est = %i et lambda = %i (%i)\n",deltat,deltat_est*alpha,lambda,deltat_est*alpha+lambda);
%            fprintf("f = %i\nt = %i\n",deltaf*Ts - (deltaf_est-phi),deltat-(deltat_est*alpha+lambda))
%             fprintf("deltat = %i est = %i \n",deltat,deltat2);
            
            
            pream_est = DemodChirp(y4,SF(k));
             
             
            y4 = Dechirp(y3(length(pream):end),SF(k),B,alpha,1);
            dnums_est = DemodChirp(y4,SF(k));
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
        TEB(i,k) = error_cnt/bit_cnt;
        TEP(i,k) = err_pck/pck_cnt;
        fprintf("TES = %i à SNR = %idB et SF = %i\n",TEB(i,k),eb_n0_dB(i),SF(k))
        if TEB(i,k) == 0
            break
        end
    end
end

%% Affichage
figure
semilogy(eb_n0_dB,TEB(:,1))
hold
for i=2:length(SF)
    semilogy(eb_n0_dB,TEB(:,i))
end
legend("SF = 9","SF = 8","SF = 7")
grid ON
ylabel("BER")
xlabel("Eb/N0 (dB)")

%% PER
figure
semilogy(eb_n0_dB,TEP(:,1))
hold
for i=2:length(SF)
    semilogy(eb_n0_dB,TEP(:,i))
end
legend("SF = 9","SF = 8","SF = 7")
grid ON
ylabel("PER")
xlabel("Eb/N0 (dB)")


