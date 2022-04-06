% Calculates differential of grad E given differential dTau.
function dwGradE = getPsiTauDeriv(W, dTau)
global gbl_f;
U = W'*O(W);
Uinv = inv(U);
Usqrtinv = sqrtm(Uinv);
tmp = dHtau(W*Uinv, dTau);
dwGradE = gbl_f*(tmp-O(W*Uinv*(W'*tmp)));
end