function merge_l_b_inputs(lake1,lake2,lake3_filename)
% function merge_l_b_inputs(lake1,lake2,lake3_filename)
% Developed by ABL, 3.11.2015.

% This MATLAB® function merge catchment input from two lakes into a new lake with calculated flow and
% concentrations based on the flow and concentrations of the two upstream lakes. 
% Lake1 input is assumed from a txt file.
% Lake2 input is assumed from a loaded mat-file. 
% Changes can be made to the function to handle several formats input data format.
% Lake3_filename is assumed a given file with .txt extensions. Adjustmest can
% be done to add the extension to the file if needed. 

disp('Merging script ...')
% opening file lake1
fid=fopen(lake1);
fgetl(fid);
tline=fgetl(fid);
%disp(tline)
header=textscan(tline,'%s','delimiter','\t');
n=length(header{1});
form=repmat('%f',1,n);
data=textscan(fid,form,'delimiter','\t');

% creating sturcture names from header 
for i=1:n
     name=header{1}{i};
     data1.(name)=data{i};
end
fclose(fid); 

k=1;
for i=11:n
    name=header{1}{i};
    data2.(name)=lake2(:,k);
    k=k+1;   
end

fac=1./(data1.InflowQ+data2.InflowQ);
nylake(1:length(data1.InflowQ),1:n)=NaN;
for i=1:10
     name=header{1}{i};
     nylake(:,i)=data{i};
end
nylake(:,11)=data1.InflowQ+data2.InflowQ;
nylake(:,12)=fac.*(data1.InflowQ.*data1.InflowT + data2.InflowQ.*data2.InflowT);
nylake(:,13)=fac.*(data1.InflowQ.*data1.InflowC + data2.InflowQ.*data2.InflowC);
nylake(:,14)=fac.*(data1.InflowQ.*data1.InflowSS + data2.InflowQ.*data2.InflowSS);
nylake(:,15)=fac.*(data1.InflowQ.*data1.InflowTP + data2.InflowQ.*data2.InflowTP);
nylake(:,16)=fac.*(data1.InflowQ.*data1.InflowDIP + data2.InflowQ.*data2.InflowDIP);
nylake(:,17)=fac.*(data1.InflowQ.*data1.InflowChla + data2.InflowQ.*data2.InflowChla);
nylake(:,18)=fac.*(data1.InflowQ.*data1.DIC + data2.InflowQ.*data2.DIC);
nylake(:,19)=fac.*(data1.InflowQ.*data1.DOC + data2.InflowQ.*data2.DOC);
nylake(:,20)=fac.*(data1.InflowQ.*data1.DO + data2.InflowQ.*data2.DO);

clear fid; 
strform=repmat('%s\t',1,n-1);
numform=repmat('%5.2f\t',1,n-1);
navn=char(header{1}{:});
fid=fopen(lake3_filename,'w');
fprintf(fid,'Mysterious useless line \n');
fprintf(fid,[strform,'%s\n'],header{1}{:});
for j=1:length(nylake)
    fprintf(fid,[numform,'%5.2f\n'],nylake(j,:));
end
fclose(fid);
end
