function [Psi, epsilon] = getPsi(W)

global gbl_f;
global gbl_kpoints;
Y = zeros(size(W));

for k = [1:gbl_kpoints]
    Wk = W(:,:,k);
    Uk = Wk'*O(Wk);
    uinv = inv(Uk)';
    usqrtinv = sqrtm(uinv);
    Y(:,:,k) = Wk*usqrtinv;
end

n = getn(Y, gbl_f);
Psi=zeros(size(Y));

for k = [1:gbl_kpoints]
    mu = Y(:,:,k)'*H(Y(:,:,k), n, k);
    [D, epsilon]=eig(mu);
    epsilon=real(diag(epsilon));
	Psi(:,:,k) = Y(:,:,k)*D;
end