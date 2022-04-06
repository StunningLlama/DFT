function grad=getgrad(W)
global gbl_f;
U = W'*O(W);
uinv = inv(U)';
usqrtinv = sqrtm(uinv);
F = diag(gbl_f);
Htilde = usqrtinv*(W'*H(W))*usqrtinv;
grad = (H(W) - O(W)*uinv*(W'*H(W)))*(usqrtinv*F*usqrtinv) + O(W)*(usqrtinv*Q(Htilde*F - F*Htilde, U));
end