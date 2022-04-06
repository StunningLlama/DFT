% Calculates differential of grad E given differential dW.
function dwGradE = getPsiPsiDeriv(W, dW)
global gbl_f;
U = W'*O(W);
Uinv = inv(U);
Usqrtinv = sqrtm(Uinv);
%F = diag(gbl_f);
%Htilde = Usqrtinv*(W'*H(W))*Usqrtinv;
dUinv= -Uinv*(dW'*O(W) + W'*O(dW))*Uinv;
tmp1 = H(W*Uinv);
% disp2(tmp1)
tmp2 = dH(W, dW*Uinv)+H2(W, dW*Uinv)+H2(W, W*dUinv); %TODO THINK
% disp("testing");
% disp2(dH(W*Uinv, dW));
% disp2(H2(W, dW*Uinv));
% disp2(H2(W, W*dUinv));
% disp2(tmp2)
dwGradE = gbl_f*((tmp2-O(W*Uinv*(W'*tmp2)) - O(dW*Uinv*(W'*tmp1) + W*dUinv*(W'*tmp1) + W*Uinv*(dW'*tmp1))));
% disp2(tmp2-O(W*Uinv*(W'*tmp2)));
% disp2( O(dW*Uinv*(W'*tmp1) + W*dUinv*(W'*tmp1) + W*Uinv*(dW'*tmp1)));

end