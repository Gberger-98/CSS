function [nums] = DeDiffp(Dnums,SF)
    nums = zeros(size(Dnums));
    for i=length(Dnums):-1:2
       nums(i)= mod(Dnums(i)-Dnums(i-1),2^SF);
    end
    nums(1)=Dnums(1);
    nums = round(nums);
end

