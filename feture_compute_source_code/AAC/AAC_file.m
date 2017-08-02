%Amino Acid Composition Descriptor
%This function calculates the Amino Acid Composition descriptor (Dim 20) by
%invoking function AAC
%The input data is a cell contains all the sequences
%Return a cell includes converted vectors
%Written by KeMeng
%Supervised by ZhangWen

function [result]=AAC_file(data)
 
 [l,m]=size(data);
 result=cell(l,1);

 for j=1:l
     str=data{j};
     result{j}=AAC(str);
 end

end