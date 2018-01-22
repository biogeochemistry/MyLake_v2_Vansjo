function res = nrmsd_mean(s, o)
    % NRMSE - normalized to the mean observed value
    res = sqrt(mean((s - o).^2))/mean(o);
end