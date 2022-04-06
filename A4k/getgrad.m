function grad=getgrad(W)
global gbl_f;
global gbl_kpoints;
global gbl_weights;
grad = zeros(size(W));
Y = zeros(size(W));
U = zeros(size(W, 2), size(W, 2), gbl_kpoints);
uinv = zeros(size(W, 2), size(W, 2), gbl_kpoints);
usqrtinv = zeros(size(W, 2), size(W, 2), gbl_kpoints);


for k = [1:gbl_kpoints]
    Wk = W(:,:,k);
    U(:,:,k) = Wk'*O(Wk);
    uinv(:,:,k) = inv(U(:,:,k))';
    usqrtinv(:,:,k) = sqrtm(uinv(:,:,k));
    Y(:,:,k) = Wk*usqrtinv(:,:,k);
end

n = getn(Y, gbl_f);

for k = [1:gbl_kpoints]
    Yk = Y(:,:,k);
    Wk = W(:,:,k);
    F = diag(gbl_f);
    Htilde = usqrtinv(:,:,k)*(Wk'*H(Wk, n, k))*usqrtinv(:,:,k);
    grad(:,:,k) = gbl_weights(k)*(...
        (H(Wk, n, k) - O(Wk)*uinv(:,:,k)*(Wk'*H(Wk, n, k)))*(usqrtinv(:,:,k)*F*usqrtinv(:,:,k)) + O(Wk)*(usqrtinv(:,:,k)*Q(Htilde*F - F*Htilde, U(:,:,k))));
end
end