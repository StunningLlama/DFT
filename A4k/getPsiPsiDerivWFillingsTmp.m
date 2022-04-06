% Calculates differential of grad E given differential dW.
function dwGradE = getPsiPsiDerivWFillingsTmp(W, dW)
global gbl_f;
F = diag(gbl_f);

U = W'*O(W);
Uinv = inv(U);
Usqrtinv = sqrtm(Uinv);
dUinv= -Uinv*(dW'*O(W) + W'*O(dW))*Uinv;
dU = dW'*O(W) + W'*O(dW);
dUsqrtinv = Q(dUinv, Uinv);
Y=W*Usqrtinv;
dY = dW*Usqrtinv+W*dUsqrtinv;

Htilde = Usqrtinv*(W'*H(W))*Usqrtinv;
dHtilde = dUsqrtinv*(W'*H(W))*Usqrtinv+Usqrtinv*(dW'*H(W))*Usqrtinv+Usqrtinv*(W'*dH(W, dW, Y, dY))*Usqrtinv+Usqrtinv*(W'*H2(W, dW*Usqrtinv+W*dUsqrtinv));

dHtilde2 = dUsqrtinv*(W'*H(W))*Usqrtinv+Usqrtinv*(dW'*H(W))*Usqrtinv;
dHtilde2 = dHtilde2+dHtilde2'+Usqrtinv*(W'*dH(W, dW, Y, dY))*Usqrtinv;

disp2(dHtilde);
disp2(dHtilde2);

tmp1 = H(W*Usqrtinv*F*Usqrtinv);
tmp2 = dH(W, dW, Y, dY)*Usqrtinv*F*Usqrtinv + H2(W*Usqrtinv*F*Usqrtinv, dW*Usqrtinv*F*Usqrtinv + W*dUsqrtinv*F*Usqrtinv+W*Usqrtinv*F*dUsqrtinv);
% disp2(dH(W*Usqrtinv*F*Usqrtinv, dW*Usqrtinv*F*Usqrtinv));
% disp2(H2(W*Usqrtinv*F*Usqrtinv, dW*Usqrtinv*F*Usqrtinv + W*dUsqrtinv*F*Usqrtinv+W*Usqrtinv*F*dUsqrtinv));
% disp2(tmp1);
% disp2(tmp2);
dwGradE = -O(dW*Uinv*(W'*tmp1) + W*dUinv*(W'*tmp1) + W*Uinv*(dW'*tmp1)) + tmp2 - O(W*Uinv*W'*tmp2) ...
    + O(dW*Usqrtinv*Q(Htilde*F-F*Htilde, U))+O(W*dUsqrtinv*Q(Htilde*F-F*Htilde, U))+O(W*Usqrtinv*Q(dHtilde*F-F*dHtilde, U))+O(W*Usqrtinv*dQ(Htilde*F-F*Htilde, U, dU));
end