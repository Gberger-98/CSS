function [DnumsM] = Diffp(numsM,SF)

    DnumsM = zeros(length(numsM)+1,1);
    for i=1:length(numsM)
        DnumsM(i+1) =mod(DnumsM(i)+numsM(i),2^SF);
    end
    DnumsM=DnumsM(2:end);
end

