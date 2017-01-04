figure(1)

% reformmating into 3D array

% big_result(1...) = TP
% big_result(2..) = Chl

%mid_run = (datenum(m_stop)-datenum(m_start))/2;
post_big_run = (squeeze(big_result(2,:,400:end)));  % maximum of , remove singleton, take the chl, all runs, remove 1 yr spin up
post_big_run(imag(post_big_run) ~= 0) = NaN; % removing imaginary numbers
post_big_run(post_big_run<0) = NaN;

post_big_vector = max(post_big_run,[],2);
post_big_vector(post_big_vector<0) = 0;

% taking the middle TP_red course (constant) out of each TP_red 20
% pack
 
% taking only the TP at constant Tred 
%TvsP = reshape(post_big_run(10:20:end),20,20); % starting from col 10, it takes every 20th value, then reshapes the 1x400 into a 20x20
TvsP = reshape(post_big_vector,20,20); % starting from col 10, it takes every 20th value, then reshapes the 1x400 into a 20x20
  
contourf(TvsP)
% chl at the end (8000)

% plot RS
response_surface_plot_vs_time
