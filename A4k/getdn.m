function dn=getdn(psi, dpsi, f)
dn = real(sum((ones(size(psi, 1), 1)*f').*(conj(dpsi).*psi+conj(psi).*dpsi), 2));
end