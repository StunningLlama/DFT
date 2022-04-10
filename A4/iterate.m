function [W,E] = iterate(iters)
global gbl_active;
%# Finite difference test
global gbl_Ns;
Ns=gbl_Ns; %# Number of states
W=(randn(length(gbl_active),Ns)+i*randn(length(gbl_active),Ns));
%more off; %# View output as it is computed
%fdtest(W);

format long
%# Converge
W=sd(W,iters); %# 20 iterations of simple sd() to get nearer to the minimum
W=W*inv(sqrtm(W'*O(W))); %# Restart as orthonormal functions
W=pccg(W,50,1); %# 50 iterations of pclm from same W
E = getE(W);
end