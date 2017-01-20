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
lake1
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
     data2.(name)=data(:,i)
end
fclose(fid); 

k=1;
for i=11:n
    name=header{1}{i};
    data2.(name)=lake2(:,k);
    k=k+1;   
end

sliced = size(data2.InflowQ, 1)-1;

fac=1./(data1.InflowQ(end - sliced:end)+data2.InflowQ);

nylake(1:length(data2.InflowQ),1:n)=NaN;
for i=1:10
     name=header{1}{i};
     nylake(:,i)=data{i};
end

nylake(:,11)=data1.InflowQ(end - sliced:end)+data2.InflowQ;
nylake(:,12)=fac.*(data1.InflowQ(end - sliced:end).*data1.InflowT(end - sliced:end) + data2.InflowQ.*data2.InflowT);
nylake(:,13)=fac.*(data1.InflowQ(end - sliced:end).*data1.InflowC(end - sliced:end) + data2.InflowQ.*data2.InflowC);
nylake(:,14)=fac.*(data1.InflowQ(end - sliced:end).*data1.InflowSS(end - sliced:end) + data2.InflowQ.*data2.InflowSS);
nylake(:,15)=fac.*(data1.InflowQ(end - sliced:end).*data1.InflowTP(end - sliced:end) + data2.InflowQ.*data2.InflowTP);
nylake(:,16)=fac.*(data1.InflowQ(end - sliced:end).*data1.InflowDIP(end - sliced:end) + data2.InflowQ.*data2.InflowDIP);
nylake(:,17)=fac.*(data1.InflowQ(end - sliced:end).*data1.InflowChla(end - sliced:end) + data2.InflowQ.*data2.InflowChla);
nylake(:,18)=fac.*(data1.InflowQ(end - sliced:end).*data1.DIC(end - sliced:end) + data2.InflowQ.*data2.DIC);
nylake(:,19)=fac.*(data1.InflowQ(end - sliced:end).*data1.DOC(end - sliced:end) + data2.InflowQ.*data2.DOC);
nylake(:,20)=fac.*(data1.InflowQ(end - sliced:end).*data1.DO(end - sliced:end) + data2.InflowQ.*data2.DO);
nylake(:,21)=fac.*(data1.InflowQ(end - sliced:end).*data1.NO3(end - sliced:end) + data2.InflowQ.*data2.NO3);
nylake(:,22)=fac.*(data1.InflowQ(end - sliced:end).*data1.NH4(end - sliced:end) + data2.InflowQ.*data2.NH4);
nylake(:,23)=fac.*(data1.InflowQ(end - sliced:end).*data1.SO4(end - sliced:end) + data2.InflowQ.*data2.SO4);
nylake(:,24)=fac.*(data1.InflowQ(end - sliced:end).*data1.Fe2(end - sliced:end) + data2.InflowQ.*data2.Fe2);
nylake(:,25)=fac.*(data1.InflowQ(end - sliced:end).*data1.Ca2(end - sliced:end) + data2.InflowQ.*data2.Ca2);
nylake(:,26)=fac.*(data1.InflowQ(end - sliced:end).*data1.pH(end - sliced:end) + data2.InflowQ.*data2.pH);
nylake(:,27)=fac.*(data1.InflowQ(end - sliced:end).*data1.CH4(end - sliced:end) + data2.InflowQ.*data2.CH4);
nylake(:,28)=fac.*(data1.InflowQ(end - sliced:end).*data1.Fe3(end - sliced:end) + data2.InflowQ.*data2.Fe3);
nylake(:,29)=fac.*(data1.InflowQ(end - sliced:end).*data1.Al3(end - sliced:end) + data2.InflowQ.*data2.Al3);
nylake(:,30)=fac.*(data1.InflowQ(end - sliced:end).*data1.SiO4(end - sliced:end) + data2.InflowQ.*data2.SiO4);
nylake(:,31)=fac.*(data1.InflowQ(end - sliced:end).*data1.SiO2(end - sliced:end) + data2.InflowQ.*data2.SiO2);
nylake(:,32)=fac.*(data1.InflowQ(end - sliced:end).*data1.diatom(end - sliced:end) + data2.InflowQ.*data2.diatom);


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
