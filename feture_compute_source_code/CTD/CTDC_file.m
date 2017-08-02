% CTD Descriptors - Composition
% This function calculates the Composition descriptor (Dim 21) by invoking
% function CTDC
% The input data is a cell contains all the sequences
% Return a cell includes converted vectors
% Written by KeMeng
% Supervised by ZhangWen

function [result]=CTDC_file(data)
 
 l=length(data);
 result=cell(l,1);

 for j=1:l
     str=data{j};
     result{j}=CTDC(str);
     
 end

end