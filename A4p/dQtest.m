%function out=dQtest(in,U, dU)
U = randn(4);
U = U*U';
dU = randn(4);
dU = 0.001*(dU*dU');
in = randn(4)+i*randn(4);
%Q(in, B+dB) - Q(in, B)
%dQ(in, B, dB)


[W,mu]=eig(U); mu=diag(mu);

[W2,mu2]=eig(U+dU); mu2=diag(mu2);


mu_n=mu*ones(1,length(mu))
mu_m=mu_n'
M=1./(sqrt(mu_n)+sqrt(mu_m))
N = 1./(mu_m-mu_n)
N(1:1+size(N,1):end) = 0


mu_n2=mu2*ones(1,length(mu2))
mu_m2=mu_n2'
M2=1./(sqrt(mu_n2)+sqrt(mu_m2))
N2 = 1./(mu_m2-mu_n2)
N2(1:1+size(N2,1):end) = 0

dmu = diag(W'*dU*W)
dmu_n = dmu*ones(1,length(mu))
dmu_m = dmu_n'

dW = W*((W'*dU*W).*N)
dM = (dmu_n./sqrt(mu_n)+dmu_m./sqrt(mu_m))./(2*(sqrt(mu_n)+sqrt(mu_m)).^2)


dmu2 = mu2-mu;
dmu_n2 = mu_n2-mu_n;

dW2 = W2-W;
dM2 = M2-M;


out= dW*((W'*in*W).*M)*W' + W*((dW'*in*W).*M)*W' + W*((W'*in*dW).*M)*W' + W*((W'*in*W).*dM)*W' + W*((W'*in*W).*M)*dW';
%end