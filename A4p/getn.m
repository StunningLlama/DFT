function n=getn(Y,f)
global gbl_weights;
global gbl_G2;
global gbl_kpoints;
global gbl_Ns;
n = zeros(size(gbl_G2));
for k = [1:gbl_kpoints]
    for j = [1:gbl_Ns]
        psikj = cI(Y{k}(:,j), k);
        n = n + gbl_weights(k)*real(f(j)*conj(psikj).*psikj);
    end
end
end