% Calculates differential of grad E given differential dTau.
function dwGradE = getPsiTauDerivWFillings(W, dTau)

global gbl_f;
global gbl_U;
global gbl_Uinv;
global gbl_Usqrtinv;
global gbl_Htilde;
global gbl_Y;
global gbl_tmp1;
global gbl_UsiFUsi;
global gbl_WUsiFUsi;
global gbl_HFmFH;
global gbl_QHFmFH;

F = diag(gbl_f);

U = gbl_U;
Uinv = gbl_Uinv;
Usqrtinv = gbl_Usqrtinv;
Htilde = gbl_Htilde;
Y = gbl_Y;
tmp1 = gbl_tmp1;
UsiFUsi = gbl_UsiFUsi;
WUsiFUsi = gbl_WUsiFUsi;
HFmFH = gbl_HFmFH;
QHFmFH = gbl_QHFmFH;

dHtilde = Usqrtinv*W'*dHtau(W*Usqrtinv, dTau);

tmp = dHtau(W*Usqrtinv*F*Usqrtinv, dTau);
dwGradE = tmp-O(W*Uinv*(W'*tmp)) + O(Y*Q(dHtilde*F-F*dHtilde, U));
end