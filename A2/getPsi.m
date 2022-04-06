function [Psi, epsilon] = getPsi(W)
Y = W*sqrtm(inv(W'*O(W)));
mu = Y'*H(Y);
[D, epsilon]=eig(mu);
epsilon=real(diag(epsilon));
Psi = Y*D;
end