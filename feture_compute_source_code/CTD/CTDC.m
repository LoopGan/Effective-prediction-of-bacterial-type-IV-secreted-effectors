% CTD Descriptors - Composition
%This function calculates the Composition descriptor (Dim 21)
%The input parameter str represents a string, for example: 'SAFLGIPFAEPPMGPRRFLPPEPKQ'
%Return : a length 21 vector
%Written by KeMeng
%Supervised by ZhangWen

function [ctdc]=CTDC(str)
a1=0;a2=0;a3=0;
len=length(str);
result=zeros(1,21);
group.hydrophobicity='RKEDQN';
group.normwaalsvolume='GASTPDC';
group.polarity='LIFWCMVY';
group.polarizability='GASDT';
group.charge= 'KR';
group.secondarystruct='EALMQKRH';
group.solventaccess='ALFCGIVW';
  
group(2).hydrophobicity = 'GASTPHY';
group(2).normwaalsvolume='NVEQIL';
group(2).polarity='PATGS';
group(2).polarizability='CPNVEQIL';
group(2).charge='ANCQGHILMFPSTWYV';
group(2).secondarystruct='VIYCWFT';
group(2).solventaccess='RKQEND';
  
group(3).hydrophobicity = 'CLVIMFW';
group(3).normwaalsvolume='MHKFRYW';
group(3).polarity='HQRKNED';
group(3).polarizability='KMHFRYW';
group(3).charge='DE';
group(3).secondarystruct='GNPSD';
group(3).solventaccess= 'MSPTHY';

for i=1:7
switch i 
    case 1 
        tt=str;
        for t=1:len
            if strfind(group(1).hydrophobicity,tt(t))~=0
                a1=a1+1;
            else if strfind(group(2).hydrophobicity,tt(t))~=0
                a2=a2+1;
                else if strfind(group(3).hydrophobicity,tt(t))~=0 
                        a3=a3+1;
                    end
                end
            end
        end
        result(i)=a1/len;a1=0;
        result(i+1)=a2/len;a2=0;
        result(i+2)=a3/len;a3=0;
     case 2  
         tt=str;
         for t=1:len
            if strfind(group(1).normwaalsvolume,tt(t))~=0
                a1=a1+1;
            else if strfind(group(2).normwaalsvolume,tt(t))~=0
                a2=a2+1;
                else if strfind(group(3).normwaalsvolume,tt(t))~=0
                        a3=a3+1;
                    end
                end
            end
        end
        result(4)=a1/len;a1=0;
        result(5)=a2/len;a2=0;
        result(6)=a3/len;a3=0;
     case 3  
         tt=str;
         for t=1:len
            if strfind(group(1).polarity,tt(t))~=0
                a1=a1+1;
            else if strfind(group(2).polarity,tt(t))~=0
                a2=a2+1;
                else if strfind(group(3).polarity,tt(t))~=0
                      a3=a3+1;
                    end
                end
            end
        end
        result(7)=a1/len;a1=0;
        result(8)=a2/len;a2=0;
        result(9)=a3/len;a3=0;
     case 4 
         tt=str;
         for t=1:len
            if strfind(group(1).polarizability,tt(t))~=0
                a1=a1+1;
            else if strfind(group(2).polarizability,tt(t))~=0
                a2=a2+1;
                else if strfind(group(3).polarizability,tt(t))~=0
                         a3=a3+1;
                    end
                end
            end
        end
        result(10)=a1/len;a1=0;
        result(11)=a2/len;a2=0;
        result(12)=a3/len;a3=0;
     case 5  
         tt=str;
         for t=1:len
            if strfind(group(1).charge,tt(t))~=0
                a1=a1+1;
            else if strfind(group(2).charge,tt(t))~=0
                a2=a2+1;
                else if strfind(group(3).charge,tt(t))~=0
                        a3=a3+1;
                    end
                end
            end
        end
        result(13)=a1/len;a1=0;
        result(14)=a2/len;a2=0;
        result(15)=a3/len;a3=0;
     case 6  
         tt=str;
         for t=1:len
            if strfind(group(1).secondarystruct,tt(t))~=0
                a1=a1+1;
            else if strfind(group(2).secondarystruct,tt(t))~=0
                a2=a2+1;
                else if strfind(group(3).secondarystruct,tt(t))~=0
                        a3=a3+1;
                    end
                end
            end
        end
        result(16)=a1/len;a1=0;
        result(17)=a2/len;a2=0;
        result(18)=a3/len;a3=0;
     case 7 
         tt=str;
         for t=1:len
            if strfind(group(1).solventaccess,tt(t))~=0
                a1=a1+1;
            else if strfind(group(2).solventaccess,tt(t))~=0
                a2=a2+1;
                else if strfind(group(3).solventaccess,tt(t))~=0
                        a3=a3+1;
                    end
                end
            end
        end
        result(19)=a1/len;a1=0;
        result(20)=a2/len;a2=0;
        result(21)=a3/len;a3=0;
end
end
ctdc=result;
end
