
%% Create lake physics output from MyLake to BROM via FABM
%% day	depth(m)	temp(deg.c)	salinity(s.u.)	density	kz(m2 s-1)
% would be great to have matlab datum instead of starting dates at 1 ...

for fabm_days = 1:(days)
stamp = ones(length(Tzt(:,fabm_days)),1); 
day_stamp = fabm_days * stamp;
salinity_stamp = 0 * stamp;
density_stamp = stamp;

fabm_depth = [1:(length(Tzt(:,fabm_days)))]';

T_mod_fabm = [Tzt(:,fabm_days)]; % all depths for day 1
Turbulence = [Kzt(:,fabm_days)]; % Kz in m2 day-1

out = [day_stamp fabm_depth T_mod_fabm salinity_stamp density_stamp Turbulence];
dlmwrite('MyLake_phys_out.txt', out, 'delimiter', '\t', 'precision', '%5.2f', '-append');
 
 end  
