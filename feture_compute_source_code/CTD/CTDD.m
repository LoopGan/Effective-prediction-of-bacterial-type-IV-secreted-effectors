% CTD Descriptors - Distribution
%This function calculates the Distribution descriptor (Dim 105)
%The input parameter str represents a string, for example: 'SAFLGIPFAEPPMGPRRFLPPEPKQ'
%Return : a length 105 vector
%Written by KeMeng
%Supervised by ZhangWen

function [ctdd]=CTDD(str)
a1=0;a2=0;a3=0;
len=length(str);
ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
result=zeros(1,105);

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
        a1=0;a2=0;a3=0;
        for t=1:len
            if strfind(group(1).hydrophobicity,tt(t))~=0
                a1=a1+1;ll1(a1)=t;
            else if strfind(group(2).hydrophobicity,tt(t))~=0
                a2=a2+1;ll2(a2)=t;
                else if strfind(group(3).hydrophobicity,tt(t))~=0 
                        a3=a3+1;ll3(a3)=t;
                    end
                end
            end
        end
        
        if a1~=0
        result(15*i-14)=100*ll1(1)/len ; 
        result(15*i-13)=100*ll1(max(floor(a1*0.25),1))/len;
        result(15*i-12)=100*ll1(max(floor(a1*0.5),1))/len;result(15*i-11)=100*ll1(floor(max(1,a1*0.75)))/len ;
        result(15*i-10)=100*ll1(a1)/len;
        
        end
        
        if a2~=0
        result(15*i-9)=100*ll2(1)/len;  
        result(15*i-8)=100*ll2(max(1,floor(a2*0.25)))/len;result(15*i-7)=100*ll2(max(1,floor(a2*0.5)))/len;result(15*i-6)=100*ll2(max(1,floor(a2*0.75)))/len ;
        result(15*i-5)=100*ll2(a2)/len;
        
        end
        
        if a3~=0
        result(15*i-4)=100*ll3(1)/len; 
        result(15*i-3)=100*ll3(max(1,floor(a3*0.25)))/len;result(15*i-2)=100*ll3(max(1,floor(a3*0.5)))/len;result(15*i-1)=100*ll3(max(1,floor(a3*0.75)))/len ;
        result(15*i)=100*ll3(a3)/len;
        end
        ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
     case 2  
         tt=str;
         a1=0;a2=0;a3=0;
         for t=1:len
            if strfind(group(1).normwaalsvolume,tt(t))~=0
                a1=a1+1;ll1(a1)=t;
            else if strfind(group(2).normwaalsvolume,tt(t))~=0
                a2=a2+1;ll2(a2)=t;
                else if strfind(group(3).normwaalsvolume,tt(t))~=0 
                        a3=a3+1;ll3(a3)=t;
                    end
                end
            end
        end
         if a1~=0
        result(15*i-14)=100*ll1(1)/len ; 
        result(15*i-13)=100*ll1(max(floor(a1*0.25),1))/len;
        result(15*i-12)=100*ll1(max(floor(a1*0.5),1))/len;result(15*i-11)=100*ll1(floor(max(1,a1*0.75)))/len ;
        result(15*i-10)=100*ll1(a1)/len;
      
        end
        
        if a2~=0
        result(15*i-9)=100*ll2(1)/len;  
        result(15*i-8)=100*ll2(max(1,floor(a2*0.25)))/len;result(15*i-7)=100*ll2(max(1,floor(a2*0.5)))/len;result(15*i-6)=100*ll2(max(1,floor(a2*0.75)))/len ;
        result(15*i-5)=100*ll2(a2)/len;
        
        end
        
        if a3~=0
        result(15*i-4)=100*ll3(1)/len; 
        result(15*i-3)=100*ll3(max(1,floor(a3*0.25)))/len;result(15*i-2)=100*ll3(max(1,floor(a3*0.5)))/len;result(15*i-1)=100*ll3(max(1,floor(a3*0.75)))/len ;
        result(15*i)=100*ll3(a3)/len;
        end
        ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
     case 3  
         tt=str;
         a1=0;a2=0;a3=0;
         for t=1:len
             if strfind(group(1).polarity,tt(t))~=0
                a1=a1+1;ll1(a1)=t;
            else if strfind(group(2).polarity,tt(t))~=0
                a2=a2+1;ll2(a2)=t;
                else if strfind(group(3).polarity,tt(t))~=0 
                        a3=a3+1;ll3(a3)=t;
                    end
                end
            end
        end
        if a1~=0
        result(15*i-14)=100*ll1(1)/len ; 
        result(15*i-13)=100*ll1(max(floor(a1*0.25),1))/len;
        result(15*i-12)=100*ll1(max(floor(a1*0.5),1))/len;result(15*i-11)=100*ll1(floor(max(1,a1*0.75)))/len ;
        result(15*i-10)=100*ll1(a1)/len;
       
        end
        
        if a2~=0
        result(15*i-9)=100*ll2(1)/len;  
        result(15*i-8)=100*ll2(max(1,floor(a2*0.25)))/len;result(15*i-7)=100*ll2(max(1,floor(a2*0.5)))/len;result(15*i-6)=100*ll2(max(1,floor(a2*0.75)))/len ;
        result(15*i-5)=100*ll2(a2)/len;
        
        end
        
        if a3~=0
        result(15*i-4)=100*ll3(1)/len; 
        result(15*i-3)=100*ll3(max(1,floor(a3*0.25)))/len;result(15*i-2)=100*ll3(max(1,floor(a3*0.5)))/len;result(15*i-1)=100*ll3(max(1,floor(a3*0.75)))/len ;
        result(15*i)=100*ll3(a3)/len;
        end
        ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
     case 4 
         tt=str;a1=0;a2=0;a3=0;
         for t=1:len
            if strfind(group(1).polarizability,tt(t))~=0
                a1=a1+1;ll1(a1)=t;
            else if strfind(group(2).polarizability,tt(t))~=0
                a2=a2+1;ll2(a2)=t;
                else if strfind(group(3).polarizability,tt(t))~=0 
                        a3=a3+1;ll3(a3)=t;
                    end
                end
            end
        end
         if a1~=0
        result(15*i-14)=100*ll1(1)/len ; 
        result(15*i-13)=100*ll1(max(floor(a1*0.25),1))/len;
        result(15*i-12)=100*ll1(max(floor(a1*0.5),1))/len;result(15*i-11)=100*ll1(floor(max(1,a1*0.75)))/len ;
        result(15*i-10)=100*ll1(a1)/len;
        
        end
        
        if a2~=0
        result(15*i-9)=100*ll2(1)/len;  
        result(15*i-8)=100*ll2(max(1,floor(a2*0.25)))/len;result(15*i-7)=100*ll2(max(1,floor(a2*0.5)))/len;result(15*i-6)=100*ll2(max(1,floor(a2*0.75)))/len ;
        result(15*i-5)=100*ll2(a2)/len;
      
        end
        
        if a3~=0
        result(15*i-4)=100*ll3(1)/len; 
        result(15*i-3)=100*ll3(max(1,floor(a3*0.25)))/len;result(15*i-2)=100*ll3(max(1,floor(a3*0.5)))/len;result(15*i-1)=100*ll3(max(1,floor(a3*0.75)))/len ;
        result(15*i)=100*ll3(a3)/len;
        end
        ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
     case 5  
         tt=str;a1=0;a2=0;a3=0;
         for t=1:len
            if strfind(group(1).charge,tt(t))~=0
                a1=a1+1;ll1(a1)=t;
            else if strfind(group(2).charge,tt(t))~=0
                a2=a2+1;ll2(a2)=t;
                else if strfind(group(3).charge,tt(t))~=0 
                        a3=a3+1;ll3(a3)=t;
                    end
                end
            end
        end
       if a1~=0
        result(15*i-14)=100*ll1(1)/len ; 
        result(15*i-13)=100*ll1(max(floor(a1*0.25),1))/len;
        result(15*i-12)=100*ll1(max(floor(a1*0.5),1))/len;result(15*i-11)=100*ll1(floor(max(1,a1*0.75)))/len ;
        result(15*i-10)=100*ll1(a1)/len;
        
        end
        
        if a2~=0
        result(15*i-9)=100*ll2(1)/len;  
        result(15*i-8)=100*ll2(max(1,floor(a2*0.25)))/len;result(15*i-7)=100*ll2(max(1,floor(a2*0.5)))/len;result(15*i-6)=100*ll2(max(1,floor(a2*0.75)))/len ;
        result(15*i-5)=100*ll2(a2)/len;
       
        end
        
        if a3~=0
        result(15*i-4)=100*ll3(1)/len; 
        result(15*i-3)=100*ll3(max(1,floor(a3*0.25)))/len;result(15*i-2)=100*ll3(max(1,floor(a3*0.5)))/len;result(15*i-1)=100*ll3(max(1,floor(a3*0.75)))/len ;
        result(15*i)=100*ll3(a3)/len;
        end
        ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
     case 6  
         tt=str;a1=0;a2=0;a3=0;
         for t=1:len
            if strfind(group(1).secondarystruct,tt(t))~=0
                a1=a1+1;ll1(a1)=t;
            else if strfind(group(2).secondarystruct,tt(t))~=0
                a2=a2+1;ll2(a2)=t;
                else if strfind(group(3).secondarystruct,tt(t))~=0 
                        a3=a3+1;ll3(a3)=t;
                    end
                end
            end
        end
        if a1~=0
        result(15*i-14)=100*ll1(1)/len ; 
        result(15*i-13)=100*ll1(max(floor(a1*0.25),1))/len;
        result(15*i-12)=100*ll1(max(floor(a1*0.5),1))/len;result(15*i-11)=100*ll1(floor(max(1,a1*0.75)))/len ;
        result(15*i-10)=100*ll1(a1)/len;
        
        end
        
        if a2~=0
        result(15*i-9)=100*ll2(1)/len;  
        result(15*i-8)=100*ll2(max(1,floor(a2*0.25)))/len;result(15*i-7)=100*ll2(max(1,floor(a2*0.5)))/len;result(15*i-6)=100*ll2(max(1,floor(a2*0.75)))/len ;
        result(15*i-5)=100*ll2(a2)/len;
       
        end
        
        if a3~=0
        result(15*i-4)=100*ll3(1)/len; 
        result(15*i-3)=100*ll3(max(1,floor(a3*0.25)))/len;result(15*i-2)=100*ll3(max(1,floor(a3*0.5)))/len;result(15*i-1)=100*ll3(max(1,floor(a3*0.75)))/len ;
        result(15*i)=100*ll3(a3)/len;
        end
        ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
     case 7 
         tt=str;a1=0;a2=0;a3=0;
         for t=1:len
            if strfind(group(1).solventaccess,tt(t))~=0
                a1=a1+1;ll1(a1)=t;
            else if strfind(group(2).solventaccess,tt(t))~=0
                a2=a2+1;ll2(a2)=t;
                else if strfind(group(3).solventaccess,tt(t))~=0 
                        a3=a3+1;ll3(a3)=t;
                    end
                end
            end
        end
         if a1~=0
        result(15*i-14)=100*ll1(1)/len ; 
        result(15*i-13)=100*ll1(max(floor(a1*0.25),1))/len;
        result(15*i-12)=100*ll1(max(floor(a1*0.5),1))/len;result(15*i-11)=100*ll1(floor(max(1,a1*0.75)))/len ;
        result(15*i-10)=100*ll1(a1)/len;
        
        end
        
        if a2~=0
        result(15*i-9)=100*ll2(1)/len;  
        result(15*i-8)=100*ll2(max(1,floor(a2*0.25)))/len;result(15*i-7)=100*ll2(max(1,floor(a2*0.5)))/len;result(15*i-6)=100*ll2(max(1,floor(a2*0.75)))/len ;
        result(15*i-5)=100*ll2(a2)/len;
        
        end
        
        if a3~=0
        result(15*i-4)=100*ll3(1)/len; 
        result(15*i-3)=100*ll3(max(1,floor(a3*0.25)))/len;result(15*i-2)=100*ll3(max(1,floor(a3*0.5)))/len;result(15*i-1)=100*ll3(max(1,floor(a3*0.75)))/len ;
        result(15*i)=100*ll3(a3)/len;
        end
        ll1=zeros(1,len);ll2=zeros(1,len);ll3=zeros(1,len);
end
end
ctdd=result;
end
