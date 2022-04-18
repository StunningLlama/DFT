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
global gbl_kpoints;
global gbl_weights;

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

dHtilde = [];
dwGradE = [];

for k = [1:gbl_kpoints]
    dHtilde(:,:,k) = Usqrtinv(:,:,k)*W(:,:,k)'*dHtau(W(:,:,k)*Usqrtinv(:,:,k), dTau);
    
    tmp = dHtau(W(:,:,k)*Usqrtinv(:,:,k)*F*Usqrtinv(:,:,k), dTau);
    dwGradE(:,:,k) = tmp-O(W(:,:,k)*Uinv(:,:,k)*(W(:,:,k)'*tmp)) + O(Y(:,:,k)*Q(dHtilde(:,:,k)*F-F*dHtilde(:,:,k), U(:,:,k)));
    dwGradE(:,:,k) = dwGradE(:,:,k)*gbl_weights(k);
end
end