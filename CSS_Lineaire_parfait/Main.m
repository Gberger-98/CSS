clc;
% close all;
clear;

%% Paramètre

B = 125e3; %Largeur de bande
% SF = [12,11,10,9,8,7]; %Nombre de bits par symbole
% SF = 12;
SF = [9,8,7];

alpha = 1;

N = 1000; %Nombre de symbole



eb=1;
eb_n0_dB = -16:0; %Rapport signal su bruit
% eb_n0_dB = 100;
eb_n0    = 10 .^ (eb_n0_dB/10);         % Liste des Eb/N0
sigma2   = 1 ./eb_n0; % Variance du bruit complex


TEB = zeros(length(eb_n0),length(SF));   % Tableau des TEB (resultats)
TEP = zeros(length(eb_n0),length(SF));


for k=1:length(SF)
    for i = 1:length(eb_n0)
        error_cnt = 0;
        bit_cnt   = 0;
        pck_cnt   = 0;
        err_pck   = 0;
        while bit_cnt < 1e6
            %% Generation
            bitsM = randi([0,1],[N,SF(k)]);

            numsM = bi2de(bitsM);

            %% Modulation
            S = MakeChirp(SF(k),numsM,B,alpha);

            %% Canal
            


            wl = sqrt(sigma2(i)/2)*(randn(size(S))+ 1j*randn(size(S)));
            yc = S + wl;

            %% Demodulatio
            y = Dechirp(yc,SF(k),B,alpha,1,0);

            nums_est = DemodChirp(y,SF(k));
            bitsM_est = de2bi(nums_est);

            tmp = sum(bitsM_est ~= bitsM,"all");
            

            bit_cnt = bit_cnt + N*SF(k);

            error_cnt =  error_cnt + tmp;
            if tmp > 0
                err_pck = err_pck + 1;
            end
            pck_cnt = pck_cnt+1;


        end
        TEB(i,k) = error_cnt/bit_cnt;
        TEP(i,k) = err_pck/pck_cnt;
        if TEB(i,k) == 0
            break
        end
        fprintf("TES = %i à SNR = %idB et SF = %i\n",TEB(i,k),eb_n0_dB(i),SF(k))
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
xlabel("SNR (dB)")

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
xlabel("SNR (dB)")


