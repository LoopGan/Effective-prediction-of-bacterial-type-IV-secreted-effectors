% CTD Descriptors - Transition
% This function calculates the Transition descriptor (Dim 21) by invoking
% function CTDT
% The input data is a cell contains all the sequences
% Return a cell includes converted vectors
% Written by KeMeng
% Supervised by ZhangWen
function [result]=CTDT_file(data)
 
 l=length(data);
 result=cell(l,1);

 for j=1:l
     str=data{j};
     result{j}=CTDT(str);
     
 end

end