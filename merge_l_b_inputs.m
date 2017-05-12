function merge_l_b_inputs(lake1,lake2,lake3_filename, m_start, m_stop)
% function merge_l_b_inputs(lake1,lake2,lake3_filename)
% Developed by ABL, 3.11.2015.

% This MATLAB® function merge catchment input from two lakes into a new lake with calculated flow and
% concentrations based on the flow and concentrations of the two upstream lakes. 
% Lake1 input is assumed from a txt file.
% Lake2 input is assumed from a loaded mat-file. 
% Changes can be made to the function to handle several formats input data format.
% Lake3_filename is assumed a given file with .txt extensions. Adjustmest can
% be done to add the extension to the file if needed. 

disp('Merging catchment input and Storefjorden input...')
% opening file lake1
fid=fopen(lake1);
fgetl(fid);
tline=fgetl(fid);
%disp(tline)
header=textscan(tline,'%s','delimiter','\t');
n=length(header{1});
form=repmat('%f',1,n);
data=textscan(fid,form,'delimiter','\t');
fclose(fid); 

% To work with any size of the input. Count number of elements from the end of data2 (it is usually smaller than data 1 provided by INCA):
sliced = size(lake2(:,1), 1)-1;
nylake = NaN(sliced+1,n);


time_in_file = [data{1},data{2},data{3}];
tmet=datenum(time_in_file);
idx_start = find(tmet==datenum(m_start));
idx_stop = find(tmet==datenum(m_stop));


% first 10 columns are the same
for i=1:10
     name=header{1}{i};
     column = data{i}(idx_start:idx_stop);
     data1.(name)=column;
     data2.(name)=column;
     nylake(:,i) =column;
end

% last columns are different
for i=11:n
     name=header{1}{i};
     data1.(name)=data{i}(end-sliced:end);
     data2.(name)=lake2(end-sliced:end,i-10);
end


fac=1./(data1.InflowQ+data2.InflowQ);
nylake(:,11)=data1.InflowQ+data2.InflowQ;
nylake(:,12)=fac.*(data1.InflowQ.*data1.InflowT + data2.InflowQ.*data2.InflowT);
nylake(:,13)=fac.*(data1.InflowQ.*data1.InflowC + data2.InflowQ.*data2.InflowC);
nylake(:,14)=fac.*(data1.InflowQ.*data1.InflowSS + data2.InflowQ.*data2.InflowSS);
nylake(:,15)=fac.*(data1.InflowQ.*data1.InflowTP + data2.InflowQ.*data2.InflowTP);
nylake(:,16)=fac.*(data1.InflowQ.*data1.InflowDIP + data2.InflowQ.*data2.InflowDIP);
nylake(:,17)=fac.*(data1.InflowQ.*data1.InflowChla + data2.InflowQ.*data2.InflowChla);
nylake(:,18)=fac.*(data1.InflowQ.*data1.DIC + data2.InflowQ.*data2.DIC);
nylake(:,19)=fac.*(data1.InflowQ.*data1.DOC + data2.InflowQ.*data2.DOC);
nylake(:,20)=fac.*(data1.InflowQ.*data1.O + data2.InflowQ.*data2.O);
nylake(:,21)=fac.*(data1.InflowQ.*data1.NO3 + data2.InflowQ.*data2.NO3);
nylake(:,22)=fac.*(data1.InflowQ.*data1.NH4 + data2.InflowQ.*data2.NH4);
nylake(:,23)=fac.*(data1.InflowQ.*data1.SO4 + data2.InflowQ.*data2.SO4);
nylake(:,24)=fac.*(data1.InflowQ.*data1.Fe2 + data2.InflowQ.*data2.Fe2);
nylake(:,25)=fac.*(data1.InflowQ.*data1.Ca2 + data2.InflowQ.*data2.Ca2);
nylake(:,26)=fac.*(data1.InflowQ.*data1.pH + data2.InflowQ.*data2.pH);
nylake(:,27)=fac.*(data1.InflowQ.*data1.CH4 + data2.InflowQ.*data2.CH4);
nylake(:,28)=fac.*(data1.InflowQ.*data1.Fe3 + data2.InflowQ.*data2.Fe3);
nylake(:,29)=fac.*(data1.InflowQ.*data1.Al3 + data2.InflowQ.*data2.Al3);
nylake(:,30)=fac.*(data1.InflowQ.*data1.SiO4 + data2.InflowQ.*data2.SiO4);
nylake(:,31)=fac.*(data1.InflowQ.*data1.SiO2 + data2.InflowQ.*data2.SiO2);
nylake(:,32)=fac.*(data1.InflowQ.*data1.diatom + data2.InflowQ.*data2.diatom);
nylake(:,33)=fac.*(data1.InflowQ.*data1.POC + data2.InflowQ.*data2.POC);



clear fid; 
strform=repmat('%s\t',1,n-1);
numform=repmat('%5.2f\t',1,n-1);
navn=char(header{1}{:});
fid=fopen(lake3_filename,'w');
fprintf(fid,'Mysterious useless line \n');
fprintf(fid,[strform,'%s\n'],header{1}{:});
dlmwrite(lake3_filename, nylake,'delimiter', '\t', '-append','precision','%5.2f')
fclose(fid);
end
