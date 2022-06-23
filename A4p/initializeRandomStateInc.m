function W=initializeRandomStateInc(q)
global gbl_active;
global gbl_Ns;
global gbl_kpoints;
Ns=gbl_Ns; %# Number of states
W = {}
for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    
    W{Tk_k}=(randn(length(gbl_active{Tk}),Ns)+i*randn(length(gbl_active{Tk}),Ns));
    W{k_Tk}=(randn(length(gbl_active{k}),Ns)+i*randn(length(gbl_active{k}),Ns));
end
end