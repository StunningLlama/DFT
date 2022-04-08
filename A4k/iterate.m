function [W,E] = iterate(sdits)
global gbl_active;
%# Finite difference test
global gbl_Ns;
global gbl_kpoints;
Ns=gbl_Ns; %# Number of states
W=(randn(length(gbl_active),Ns, gbl_kpoints)+i*randn(length(gbl_active),Ns, gbl_kpoints));
%more off; %# View output as it is computed
%fdtest(W);

format long
%# Converge
W=sd(W,sdits); %# 20 iterations of simple sd() to get nearer to the minimum
W = orthonormalize(W); %# Restart as orthonormal functions
W=pccg(W,50,1); %# 50 iterations of pclm from same W
E = getE(W);
end