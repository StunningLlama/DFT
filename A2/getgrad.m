function grad=getgrad(W)
uinv = inv(W'*O(W))';
grad = (H(W) - O(W*uinv)*(W'*H(W)))*uinv;
end