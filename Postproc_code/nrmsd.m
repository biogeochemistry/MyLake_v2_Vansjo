function res = nrmsd(s, o)
    % NRMSE - normalized to standard deviation of RMSE
    res = sqrt(mean((s - o).^2))/std(o);
end