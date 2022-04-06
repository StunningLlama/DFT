% Calculates differential of grad E given differential dW.
function dwGradE = getPsiPsiDerivWFillings(W, dW)
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
global gbl_HWUsi;

F = diag(gbl_f);

U = gbl_U;
Uinv = gbl_Uinv;
Usqrtinv = gbl_Usqrtinv;
Htilde = gbl_Htilde;
Y = gbl_Y;
tmp1 = gbl_tmp1;
UsiFUsi = gbl_UsiFUsi; %U^(-1/2)FU^(-1/2)
WUsiFUsi = gbl_WUsiFUsi; %WU^(-1/2)FU^(-1/2)
HFmFH = gbl_HFmFH; %~HF - F~H
QHFmFH = gbl_QHFmFH; %Q(~HF - F~H)
HWUsi= gbl_HWUsi;

dUinv= -Uinv*(dW'*O(W) + W'*O(dW))*Uinv;
dU = dW'*O(W) + W'*O(dW);
dUsqrtinv = Q(dUinv, Uinv);
dY = dW*Usqrtinv+W*dUsqrtinv;

dHW = dH(W, dW, Y, dY);
dHtilde = (dUsqrtinv*W'+Usqrtinv*dW')*HWUsi;
dHtilde = dHtilde+dHtilde'+Usqrtinv*(W'*dHW)*Usqrtinv;
        %Todo think about this term.
tmp2 = dHW*UsiFUsi + H2(WUsiFUsi, dW*UsiFUsi + W*dUsqrtinv*F*Usqrtinv+W*Usqrtinv*F*dUsqrtinv);
% disp2(dH(W*Usqrtinv*F*Usqrtinv, dW*Usqrtinv*F*Usqrtinv));
% disp2(H2(W*Usqrtinv*F*Usqrtinv, dW*Usqrtinv*F*Usqrtinv + W*dUsqrtinv*F*Usqrtinv+W*Usqrtinv*F*dUsqrtinv));
% disp2(tmp1);
% disp2(tmp2);
dwGradE = -O(dW*Uinv*(W'*tmp1) + W*dUinv*(W'*tmp1) + W*Uinv*(dW'*tmp1)) + tmp2 - O(W*(Uinv*(W'*tmp2))) ...
    + O(dY*QHFmFH + Y*Q(dHtilde*F-F*dHtilde, U) + Y*dQ(HFmFH, U, dU));
end