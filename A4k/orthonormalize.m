function Y = orthonormalize(W)
global gbl_kpoints;
Y = {};
for k = [1:gbl_kpoints]
    Wk = W{k};
    Uk = Wk'*O(Wk);
    uinv = inv(Uk)';
    usqrtinv = sqrtm(uinv);
    Y{k} = Wk*usqrtinv;
end
end