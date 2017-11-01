function res = rmsd(s, o)
    res = sqrt(mean((s - o).^2));
end