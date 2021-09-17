clc;
% close all;
clear;

%% Paramètre


% SF = [12,11,10,9,8,7]; %Nombre de bits par symbole
SF = [7];
% SF = 7;



alpha = 10;
B = 125e3; %Largeur de bande
te = 1/(alpha*B);
N = 1000; %Nombre de symbole


eb=1;
eb_n0_dB = -16:1:-5; %Rapport signal su bruit
% eb_n0_dB = 100;
eb_n0    = 10 .^ (eb_n0_dB/10);         % Liste des Eb/N0
sigma2   = 1 ./eb_n0; % Variance du bruit complex


TEB = zeros(length(eb_n0),length(SF));   % Tableau des TEB (resultats)
TEP = zeros(length(eb_n0),length(SF));
%% défaut

maxcfo = B/4.2;


for k=1:length(SF)
    pream = MakePream(SF(k),B,alpha);
    Ts = 2^SF(k)/B;
    for i = 1:length(eb_n0)
        error_cnt = 0;
        bit_cnt   = 0;
        pck_cnt     = 0;
        err_pck_cnt = 0;
        while bit_cnt < 1e6
            %% Generation
            bitsM = randi([0,1],[N,SF(k)]);

            numsM = bi2de(bitsM);

            %% Modulation
            
            S = MakeChirp(SF(k),numsM,B,alpha);
            S = [pream,S];
            %% Canal
            

            deltat = randi([0,alpha*(2^SF(k))/4-1],1);
%             deltat = randi([-(alpha*(2^SF(k))/4-1),0],1);
            deltaf = randi(round(maxcfo));
            if deltat>=0
                S = [zeros(1,deltat),S];
            else
                S = S(-deltat:end);
            end
            wl = sqrt(sigma2(i)/2)*(randn(size(S))+ 1j*randn(size(S)));
            yc = S.*exp(1j*2*pi*deltaf*te*[0:length(S)-1])+ wl;

            %% Demodulation
           
            y = Dechirp(yc,SF(k),B,alpha,2,1);
            
            [phi,lambda] = synchro_fraq(y,SF(k),B,alpha);
            
            y = yc(lambda+1:end).*exp(-1j*2*pi*(-phi)*te/Ts*[0:length(yc(lambda+1:end))-1]);
            y = Dechirp(y,SF(k),B,alpha,1,1);
            [deltaf_est,deltat_est] = synchro_int(y,SF(k),B);
            
            y2 = yc.*exp(-1j*2*pi*(deltaf_est-phi)*te/Ts*[0:length(yc)-1]); %synchro freq
            if deltat_est<2^SF(k)/4
                y3 = y2(deltat_est*alpha+lambda+1:end); %Synchro temps
            else
                y3 = y2(lambda+1:end); %Synchro temps
            end

            %y3 = y2(deltat_est*alpha+1:end);
            y4 = Dechirp(y3,SF(k),B,alpha,1,1);
            
%              fprintf("delta f =  %i, est = %i, dif = %i \n",deltaf*Ts,deltaf_est-phi,deltaf*Ts - (deltaf_est-phi))
%             fprintf("delta t = %i, est = %i et lambda = %i (%i)\n",deltat,deltat_est*alpha,lambda,deltat_est*alpha+lambda);
%            fprintf("f = %i\nt = %i\n",deltaf*Ts - (deltaf_est-phi),deltat-(deltat_est*alpha+lambda))

            
            nums_est_ent = DemodChirp(y4,SF(k));

            pream_est = nums_est_ent(1:12); 
            nums_est = nums_est_ent(13:end);
            bitsM_est = de2bi(mod(nums_est-mode(pream_est(1:8)),2^SF(k)),SF(k));

            tmp = sum(bitsM_est ~= bitsM((end-length(bitsM_est)+1):end,:),"all");


%             if  tmp > 0
%                 fprintf("tmp = %i\n",tmp)
%             end
            
            bit_cnt = bit_cnt + N*SF(k);
            error_cnt =  error_cnt + tmp;
            pck_cnt     = pck_cnt + 1;
            if tmp >0
                err_pck_cnt = err_pck_cnt + 1;
            end


        end
        TEB(i,k) = error_cnt/bit_cnt;
        TEP(i,k) = err_pck_cnt/pck_cnt;
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


%%
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

