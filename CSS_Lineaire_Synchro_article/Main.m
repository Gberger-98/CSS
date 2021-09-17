clc;
% close all;
clear;

%% Paramètre

B = 125e3; %Largeur de bande
% SF = [12,11,10,9,8,7]; %Nombre de bits par symbole
SF = [9,8,7];
% SF = 7;


alpha = 10;
te = 1/(alpha*B);
% B_E = [-1, 0.15 0.5 0.7];
% B_E = -1;

N = 100; %Nombre de symbole



Nbpream = 10;

eb=1;
% eb_n0_dB = 100; %Rapport signal sur bruit
eb_n0_dB = -16:-5;
eb_n0    = 10 .^ (eb_n0_dB/10);         % Liste des Eb/N0
sigma2   = 1 ./eb_n0; % Variance du bruit complex


TEB = zeros(length(eb_n0),length(SF));   % Tableau des TEB (resultats)
PER = zeros(length(eb_n0),length(SF));
c = [];
maxcfo = B/4.2;

for k=1:length(SF)
    M=2^SF(k);
    Ts = 2^SF(k)/B;
    T=M/B;
    for i = 1:length(eb_n0)
        error_cnt   = 0;
        bit_cnt     = 0;
        pck_cnt     = 0;
        err_pck_cnt = 0;
        pream = MakeChirp(SF(k),zeros(Nbpream,1),B,alpha);
        up = pream(1:alpha*2^SF(k));
        down = conj(up);
        pream = [pream down down];
        while bit_cnt < 1e6
            %% Generation
            bitsM = randi([0,1],[N,SF(k)]);
            numsM = bi2de(bitsM);


            %% Modulation
            S = MakeChirp(SF(k),numsM,B,alpha);
            S = [pream,S];

            %% Canal

            deltat = randi([0,alpha*(2^SF(k))/4-1],1);
            deltaf = randi(round(maxcfo));
            S = [zeros(1,deltat),S];
            wl = sqrt(sigma2(i)/2)*(randn(size(S))+ 1j*randn(size(S)));

            yc = S.*exp(1j*2*pi*deltaf*te*[0:length(S)-1])+ wl;



            %% Demodulation
            y1 = Dechirp(yc,SF(k),B,alpha,2,1);
            [phi,~] = synchro_fraq(y1,SF(k),B,alpha);

            yt = yc.*exp(-1j*2*pi*(-phi)*te/T*[0:length(yc)-1]);
            y = Dechirp(yt,SF(k),B,alpha,1,1);
            lambda1 = firststo(y,SF(k));
            deb = floor(alpha*lambda1+1);
            if lambda1>0
                yt = yc(deb:end).*exp(-1j*2*pi*(-phi)*te/T*[0:length(yc(deb:end))-1]);
            else
                yt = yc.*exp(-1j*2*pi*(-phi)*te/T*[0:length(yc)-1]);
            end
            y = Dechirp(yt,SF(k),B,alpha,1,1);
            [deltaf_est,deltat_est,lambda] = finalSync(y,SF(k));
            if deb ==0
                deb=1;
            end
            if lambda1>0
                deb = abs(floor(alpha*(deltat_est-lambda+lambda1)))+1;
            else
                deb = abs(floor(alpha*(deltat_est-lambda)))+1;
            end
            yt = yc(deb:end).*exp(-1j*2*pi*(deltaf_est-phi)*te/T*[0:length(yc(deb:end))-1]);
            y4 = Dechirp(yt,SF(k),B,alpha,1,1);
            deltat2 = deb;

            deltaf2 = round((deltaf_est-phi)/Ts);
            nums_est_ent = DemodChirp(y4,SF(k));

            pream_est = nums_est_ent(1:12); 
            nums_est = nums_est_ent(13:end);
            offset = mode(pream_est(1:8));

            bitsM_est = de2bi(mod(nums_est-offset,2^SF(k)),SF(k));
            tmp = sum(bitsM_est ~= bitsM((end-length(bitsM_est)+1):end,:),"all");


            bit_cnt = bit_cnt + N*SF(k);

            error_cnt =  error_cnt + tmp;

            pck_cnt     = pck_cnt + 1;
            if tmp >0
                err_pck_cnt = err_pck_cnt + 1;
            end

        end
        TEB(i,k) = error_cnt/bit_cnt;
        PER(i,k) = err_pck_cnt/pck_cnt;

        fprintf("TEB = %1.2e à SNR = %2.1fdB et SF = %i\n",TEB(i,k),eb_n0_dB(i),SF(k))

        if TEB(i,k) == 0
            break
        end
    end
end

%% Affichage
figure
semilogy(eb_n0_dB,TEB(:,1))
hold
styl = ["--","--s","--*","--x","--d","--p","--^","-->","--<","--h"];
for i=2:length(SF)
    semilogy(eb_n0_dB,TEB(:,i))
end

%% Legend
Legend = strings(length(SF),1);

for i=1:length(SF)
    Legend(i) = sprintf("SF = %i ",SF(i)); 
end
legend(Legend);
grid ON
ylabel("BER")
xlabel("SNR (dB)")

%% PER
figure
semilogy(eb_n0_dB,PER(:,1))
hold
styl = ["--","--s","--*","--x","--d","--p","--^","-->","--<","--h"];
for i=2:length(SF)
    semilogy(eb_n0_dB,PER(:,i))
end

%% Legend
Legend = strings(length(SF),1);

for i=1:length(SF)
   Legend(i) = sprintf("SF = %i ",SF(i)); 
end
legend(Legend);
grid ON
ylabel("PER")
xlabel("SNR (dB)")

