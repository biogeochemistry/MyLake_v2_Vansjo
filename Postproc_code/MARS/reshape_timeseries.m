%% This is to reshape the entire simulation of 41 years so that each column is a year%
% input is a C&P from excel file DOC Sensitivity analysis

all_shaved=all_brut(1:14965,:)

DOC_cst_1m=all_brut(1:14965,1)
DOC_cst_8m=all_brut(1:14965,2)
Temp_cst_1m=all_brut(1:14965,3)
Temp_cst_8m=all_brut(1:14965,4)
all_cst_1m=all_brut(1:14965,5)
all_cst_8m=all_brut(1:14965,6)
%hist_cst_1m=all_brut(1:14965,7) % not for diff time series 
%hist_cst_8m=all_brut(1:14965,8) % not for diff time series 

DOC_1m_reshape=reshape(DOC_cst_1m,365,[])
DOC_8m_reshape=reshape(DOC_cst_8m,365,[])
Temp_8m_reshape=reshape(Temp_cst_8m,365,[])
Temp_1m_reshape=reshape(Temp_cst_1m,365,[])
all_cst_1m_reshape=reshape(all_cst_1m,365,[])
all_cst_8m_reshape=reshape(all_cst_8m,365,[])
%hist_cst_1m_reshape=reshape(hist_cst_1m,365,[]) % not for diff time series 
%hist_cst_8m_reshape=reshape(hist_cst_8m,365,[]) % not for diff time series 

all_daily_means= zeros(365,1)

all_daily_means(:,1)=mean(DOC_1m_reshape,2)
all_daily_means(:,2)=mean(DOC_8m_reshape,2)
all_daily_means(:,3)=mean(Temp_1m_reshape,2)
all_daily_means(:,4)=mean(Temp_8m_reshape,2)
all_daily_means(:,5)=mean(all_cst_1m_reshape,2)
all_daily_means(:,6)=mean(all_cst_8m_reshape,2)
%all_daily_means(:,7)=mean(hist_cst_1m_reshape,2)  % not for diff time series 
%all_daily_means(:,8)=mean(hist_cst_8m_reshape,2)  % not for diff time series 