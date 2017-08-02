% CTD Descriptors - Distribution
% This function calculates the Distribution descriptor (Dim 105) by
% invoking function CTDD
% The input data is a cell contains all the sequences
% Return a cell includes converted vectors
% Written by KeMeng
% Supervised by ZhangWen

function [result]=CTDD_file(data)
 
 l=length(data);
 result=cell(l,1);

 for j=1:l
     str=data{j};
     result{j}=CTDD(str);  
 end

end