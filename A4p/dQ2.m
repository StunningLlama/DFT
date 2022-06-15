function out=dQ2(in,U, dU)
[W dmu]=eig(dU);
dmu = diag(dmu);

mu = ones(size(U,1),1);
mu_n=mu*ones(1,length(mu));
mu_m=mu_n';
%W'*dU*W
%dmu = diag(W'*dU*W);
dmu_n = dmu*ones(1,length(mu));
dmu_m = dmu_n';

dM = -(dmu_n./sqrt(mu_n)+dmu_m./sqrt(mu_m))./(2*(sqrt(mu_n)+sqrt(mu_m)).^2);
W'*in*W
out= W*((W'*in*W).*dM)*W';
end