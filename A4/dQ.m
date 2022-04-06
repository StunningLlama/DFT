function out=dQ(in,U, dU)
[W,mu]=eig(U); mu=diag(mu);
mu_n=mu*ones(1,length(mu));
mu_m=mu_n';
M=1./(sqrt(mu_n)+sqrt(mu_m));
N = 1./(mu_m-mu_n);
N(1:1+size(N,1):end) = 0;

dmu = diag(W'*dU*W);
dmu_n = dmu*ones(1,length(mu));
dmu_m = dmu_n';

dW = W*((W'*dU*W).*N);
dM = -(dmu_n./sqrt(mu_n)+dmu_m./sqrt(mu_m))./(2*(sqrt(mu_n)+sqrt(mu_m)).^2);
out= dW*((W'*in*W).*M)*W' + W*((dW'*in*W).*M)*W' + W*((W'*in*dW).*M)*W' + W*((W'*in*W).*dM)*W' + W*((W'*in*W).*M)*dW';
end