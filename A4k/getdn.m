function dn=getdn(psi, dpsi, f)
global gbl_G2;
dn = zeros(size(gbl_G2));
global gbl_kpoints;
global gbl_weights;
for k = [1:gbl_kpoints]
    dn = dn + gbl_weights(k)*real(sum((ones(size(psi{k}, 1), 1)*f').*(conj(dpsi{k}).*psi{k}+conj(psi{k}).*dpsi{k}), 2));
end
end