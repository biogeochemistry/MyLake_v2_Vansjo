fid = fopen('IO/mylake_initial_concentrations.txt','wt');
fprintf(fid, '\n Z (m)\t  Az (m2)\t Tz (deg C)\t  Cz(kg/m3)\t  Sz (kg/m3)\t  TPz (mg/m3)\t DOPz (mg/m3)\t    Chlaz (mg/m3)\t   DOCz (mg/m3)\t    TPz_sed (mg/m3)\t Chlaz_sed (mg/m3)\t   Fvol_IM (m3/m3, dry w.)\t Hice (m)\t Hsnow (m)\t O2z (mg/m3)\t DICz (mg/m3)\t NO3z (mg/m3)\t NH4z (mg/m3)\t SO4z (mg/m3)\t HSz (mg/m3)\t H2Sz (mg/m3)\t Fe2z (mg/m3)\t Ca2z (mg/m3)\t pHz (mg/m3)\t CH4z (mg/m3)\t Fe3z (mg/m3)\t Al3z (mg/m3)\t SiO4z (mg/m3)\t SiO2z (mg/m3)\t diatomz (mg/m3)\t POCz (mg/m3)\n');
fclose(fid);


dlmwrite('IO/mylake_initial_concentrations.txt', ...
 [MyLake_results.z(:,end), ...
 MyLake_results.params.Az(:,end), ...
 MyLake_results.Tzt(:,end), ...
 MyLake_results.Czt(:,end), ...
 MyLake_results.Szt(:,end), ...
 MyLake_results.Pzt(:,end)+MyLake_results.PPzt(:,end)+MyLake_results.DOPzt(:,end)+MyLake_results.Chlzt(:,end)+MyLake_results.Czt(:,end), ...
 MyLake_results.DOPzt(:,end), ...
 MyLake_results.Chlzt(:,end), ...
 MyLake_results.DOCzt(:,end), ...
 zeros(length(MyLake_results.z),1), ...  % In_TPz_sed
 zeros(length(MyLake_results.z),1), ...  % In_Chlz_sed
 zeros(length(MyLake_results.z),1), ...   % MyLake_results.FIM(:,end)
 zeros(length(MyLake_results.z),2), ...  % ice0
 MyLake_results.O2zt(:,end), ...
 MyLake_results.DICzt(:,end), ...
 MyLake_results.NO3zt(:,end), ...
 MyLake_results.NH4zt(:,end), ...
 MyLake_results.SO4zt(:,end), ...
 MyLake_results.HSzt(:,end), ...
 MyLake_results.H2Szt(:,end), ...
 MyLake_results.Fe2zt(:,end), ...
 MyLake_results.Ca2zt(:,end), ...
 MyLake_results.pHzt(:,end), ...
 MyLake_results.CH4zt(:,end), ...
 MyLake_results.Fe3zt(:,end), ...
 MyLake_results.Al3zt(:,end), ...
 MyLake_results.SiO4zt(:,end), ...
 MyLake_results.SiO2zt(:,end), ...
 MyLake_results.diatomzt(:,end), ...
 MyLake_results.POCzt(:,end)], 'delimiter', '\t', '-append')


