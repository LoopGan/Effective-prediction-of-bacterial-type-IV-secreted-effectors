function [] = featurecomput(infile,dbfile,pssmfile,AAindexfile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% count the number of input sequences
fid1=fopen(infile,'r');
seqnum=0;
while ~feof(fid1)
    lines=fgetl(fid1);
    if lines(1)=='>'
        seqnum=seqnum+1;
    end
end
fclose(fid1);

%% Calculate PSSM
pssm(infile,dbfile,seqnum);
%% Calculate PSSM composition features
pssmcomp(pssmfile,seqnum);
%% calculate AAindex
AAindex(infile,AAindexfile);

end

function [] = pssm(infile,dbfile,seqnum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
inname=infile;
dbname=dbfile;
mkdir .\seq;
cmd1=['perl seqdiv.pl ' inname ' .\seq'];
dos(cmd1);
for num=1:seqnum
    seqname=['sp_' num2str(num) '.fa'];
    outname=['seq' num2str(num) '.pssm'];
    cmd2=['.\blast\bin\blastpgp -i .\seq\',seqname, ' -d', dbname ,' -j 3 -h 0.001 -a 8 -Q .\seq\',outname];
    dos(cmd2); 
end
end

%%
function [ ] = pssmcomp(outfile,num)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
aa='ARNDCQEGHILKMFPSTWYV';
fid1=fopen(outfile,'w');
for i=1:num
    seqfile=['./seq/','sp_',num2str(i),'.fa'];
    fid2=fopen(seqfile,'r');
    seqname=fgetl(fid2);
    seq=fgetl(fid2);
    len=length(seq);
    fclose(fid2);
    %
    infile=['./seq/','seq',num2str(i),'.pssm'];
    indata=importdata(infile);
    data1=indata.data;
    [p,q]=size(data1);
    pssm=data1(1:(p-5),1:20);
    %{
    data2=indata.textdata;
    [m,n]=size(data2);
    seq=data2(4:(m-5),2);
    %}
     % start ac
    output = [];
     for j=1:20
         pos=strfind(seq,aa(j));
         posdata=pssm(pos,:);
         if length(pos)==1
             d(1:20)=0;
             posdata=[posdata;d];
                 
         end
                 
         t=sum(posdata);
        
         if length(pos)==0
             t(1:20)=0;
         end
       
         output=[output,t];
     end
   output1=output/len;
%  output=1./(1+exp(-output1));
   %
    fprintf(fid1,'%f ',output1);
    fprintf(fid1,'\n');
end
fclose(fid1);
end

%%
function AAindex=Run_hyy(fasta_filename,feature_filename)

warning off
[header,sequence]=fastaread(fasta_filename);
onebound=numel(header);
AAindex=cell(onebound,3);
only_id=RungetonlyID(header);AAindex(:,1)=only_id;

[numeric,txt,~]=xlsread(feature_filename);
txt=char(txt);txt=txt';

for ii=1:onebound
   only_seq=sequence{ii};
   only_fea=humap(only_seq,txt,numeric);
   only_mean=mean(only_fea,2);
   AAindex{ii,2}=only_mean;
end

end

function list2=humap(list,origin,newcolumn)
[list,vSort]=sort(list);
[~,vRev]=sort(vSort);

[list,vEnd]=unique(list);
vLen=vEnd-[0,vEnd(1:end-1)];

[~,~,id]=intersect(list,origin);
V=zeros(1,sum(vLen));L=1;
for I=1:numel(id)
    v=L:L+vLen(I)-1;
    V(1,v)=id(I);
    L=L+vLen(I);
end
V=V(vRev);
list2=newcolumn(:,V);
end

function onlyID=RungetonlyID(header)
onebound=numel(header);
onlyID=cell(onebound,1);
for ii=1:onebound
    str=header{ii};
    v=find(str=='|');
    onlyID{ii}=str(1,v(1)+1:v(2)-1);
end
end


