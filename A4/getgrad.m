function grad=getgrad(W)
global gbl_f;
U = W'*O(W);
uinv = inv(U)';
usqrtinv = sqrtm(uinv);
F = diag(gbl_f);
WdagHW = W'*H(W);
Htilde = usqrtinv*WdagHW*usqrtinv;
grad = (H(W) - O(W)*uinv*(WdagHW))*(usqrtinv*F*usqrtinv) + O(W)*(usqrtinv*Q(Htilde*F - F*Htilde, U));
end