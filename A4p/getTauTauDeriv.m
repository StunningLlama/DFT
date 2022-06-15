function dXdXE = getTauTauDeriv(dXa, dXb)
global gbl_f;
global gbl_G;
global gbl_Vps;
global gbl_n;
global gbl_X;
G = gbl_G;
n = gbl_n;
dSf = getdsqSf(gbl_X, dXa, dXb);
dVtilde = gbl_Vps.*dSf;
dVtilde(1)=0.;
dVdual=cJ(dVtilde);
dXdXE = real(dVdual')*n;
end