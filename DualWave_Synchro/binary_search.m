function   index_max=binary_search(buffer,index,M,Nb_preambule,alpha,max_iter,chirp_brut)

if alpha>1
    tempLoRa = buffer(1:alpha:end);
    buffer = [];
    buffer = tempLoRa;
    
    tempChirpLoRa = chirp_brut(1:alpha:end);
    chirp_brut = [];
    chirp_brut = tempChirpLoRa;
end 
  
Nc=floor(length(buffer)/(M));
chirp_mat=repmat(chirp_brut',1, Nc);
buffer_unspread = buffer(1:Nc*M).* chirp_mat(1:end); % de-chirped buffer
    if index<=0
        index=1;
    end
 
    a=M/2;
    ind_ap=index*M+a;
    ind_av=index*M-a;
    index_max=index*M;
    ii=0;
 while ii<max_iter+1
    buffer_test_ap= buffer_unspread(1+ind_ap:ind_ap+(Nb_preambule+2)*M);
    buffer_test_av= buffer_unspread(1+ind_av:ind_av+(Nb_preambule+2)*M);
    
%     buffer_test_ap= [buffer_unspread(1+ind_ap:ind_ap+M) buffer_unspread(1+ind_ap+(Nb_preambule)*M:ind_ap+(Nb_preambule+2)*M)];
%     buffer_test_av= [buffer_unspread(1+ind_av:ind_av+M) buffer_unspread(1+ind_av+(Nb_preambule)*M:ind_av+(Nb_preambule+2)*M)];
%     
    
    S_ap = spectrogram(buffer_test_ap/sqrt(M),ones(1,M),0,2*M);
    S_av = spectrogram(buffer_test_av/sqrt(M),ones(1,M),0,2*M);
%     power_ap= max(sum(abs(S_ap),2));
%     power_av= max(sum(abs(S_av),2));
  
    power_ap= sum(max(abs(S_ap)));
    power_av= sum(max(abs(S_av)));

    if power_ap> power_av
      ind_av=index_max;
      index_max =  (ind_ap+index_max)/2;
    else
      ind_ap=index_max;  
      index_max =  (index_max+ind_av)/2;
    end
  ii=ii+1;
 end
 index_max=round(index_max);
end

