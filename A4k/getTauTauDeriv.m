function dXdXE = getTauTauDeriv(W, X, dXa, dXb)
global gbl_f;
global gbl_G;
global gbl_Vps;
G = gbl_G;
Y=W*inv(sqrtm(W'*O(W)));
n = getn(Y, gbl_f);
dSf = getdsqSf(X, dXa, dXb);
dVtilde = gbl_Vps.*dSf;
dVtilde(1)=0.;
dVdual=cJ(dVtilde);
dXdXE = real(dVdual')*n;
end