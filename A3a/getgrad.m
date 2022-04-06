function grad=getgrad(W)
global gbl_f;
uinv = inv(W'*O(W))';
grad = gbl_f*(H(W) - O(W*uinv)*(W'*H(W)))*uinv;
end