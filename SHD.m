function [output] = SHD(ch1,ch2)
%Structure Hamming Distance
%INPUT: two charateristic imsets
%OUTPUT: SHD

equal = sum(ch1==ch2);
output=length(ch1)-equal;

end