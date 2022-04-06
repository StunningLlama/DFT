function n=getn(psi,f)
n = real(sum((ones(size(psi, 1), 1)*f').*conj(psi).*psi, 2));
end