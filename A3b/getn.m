function n=getn(psi,f)
n = f*sum(conj(psi).*psi, 2);
end