function W=initializeRandomState()
global gbl_active;
global gbl_Ns;
global gbl_kpoints;
Ns=gbl_Ns; %# Number of states
W = {}
for k = [1:gbl_kpoints]
    W{k}=(randn(length(gbl_active{k}),Ns)+i*randn(length(gbl_active{k}),Ns));
end
end