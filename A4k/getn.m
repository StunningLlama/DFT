function n=getn(Y,f)
global gbl_weights;
global gbl_G2;
global gbl_kpoints;
n = zeros(size(gbl_G2));
for k = [1:gbl_kpoints]
psik = cI(Y(:,:,k));
n = n + gbl_weights(k)*real(sum((ones(size(psik, 1), 1)*f').*conj(psik).*psik, 2));
end
end