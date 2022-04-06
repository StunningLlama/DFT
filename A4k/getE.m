function E=getE(W)
global gbl_f;
global gbl_kpoints;
global gbl_Vdual;
global gbl_weights;
Y = zeros(size(W));

for k = [1:gbl_kpoints]
    Wk = W(:,:,k);
    Uk = Wk'*O(Wk);
    uinv = inv(Uk)';
    usqrtinv = sqrtm(uinv);
    Y(:,:,k) = Wk*usqrtinv;
end

n = getn(Y, gbl_f);

E = real(gbl_Vdual'*n + 0.5*n'*cJdag(O(-4*pi*Linv(O(cJ(n))))) + n'*cJdag(O(cJ(excVWN(n)))));
for k = [1:gbl_kpoints]
    E = E+gbl_weights(k)*real(-0.5*trace(diag(gbl_f)*(Y(:,:,k)'*L(Y(:,:,k), k))));
end
end