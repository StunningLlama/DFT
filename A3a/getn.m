function n=getn(psi,f)
n = sum(conj(psi).*psi, 2);
end