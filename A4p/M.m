function psiout = M(gvec, psi)
global gbl_r;
garray = ones(size(gbl_r,1),1)*gvec;
phases = exp(i*sum(gbl_r.*garray,2));
psiout = psi.*phases;
end