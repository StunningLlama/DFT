function n=getn(Y,f)
global gbl_weights;
global gbl_G2;
global gbl_kpoints;
n = zeros(size(gbl_G2));
for k = [1:gbl_kpoints]
    for j = [1:size(Y,2)]
        psikj = cI(Y(:,j,k));
        n = n + gbl_weights(k)*real(f(j)*conj(psikj).*psikj);
    end
end
end