%Amino Acid Composition Descriptor
%This function calculates the Amino Acid Composition descriptor (Dim 20)
%The input parameter str represents a string, for example: 'SAFLGIPFAEPPMGPRRFLPPEPKQ'
%Return AAC: a length 20 vector
%Written by KeMeng
%Supervised by ZhangWen

function [AAC]=AAC(str)
len=length(str);
C={'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
l=20;
result=zeros(1,l);
Result=zeros(1,l);
for i=1:len
    for t=1:l
        if strcmp(str(i) , C(t)) 
           result(t)=result(t)+1 ;
        end
    end   
end
for t=1:l
    Result(t)=result(t)/len;
end
AAC=Result;
end