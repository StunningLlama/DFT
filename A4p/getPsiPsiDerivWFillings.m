% Calculates differential of grad E given differential dW.
function dwGradE = getPsiPsiDerivWFillings(dW)
global gbl_f;
global gbl_U;
global gbl_Uinv;
global gbl_Usqrtinv;
global gbl_Htilde;
global gbl_Y;
global gbl_tmp1;
global gbl_UsiFUsi;
global gbl_WUsiFUsi;
global gbl_HFmFH;
global gbl_QHFmFH;
global gbl_HWUsi;
global gbl_IY;
global gbl_IW;
global gbl_kpoints;
global gbl_weights;
global gbl_W;
global gbl_excpn;
global gbl_excppn;
global gbl_JdagOJn;
global gbl_hsum;

F = diag(gbl_f);

U = gbl_U;
Uinv = gbl_Uinv;
Usqrtinv = gbl_Usqrtinv;
Htilde = gbl_Htilde;
Y = gbl_Y;
IW = gbl_IY;
IW = gbl_IW;
tmp1 = gbl_tmp1;
UsiFUsi = gbl_UsiFUsi; %U^(-1/2)FU^(-1/2)
WUsiFUsi = gbl_WUsiFUsi; %WU^(-1/2)FU^(-1/2)
HFmFH = gbl_HFmFH; %~HF - F~H
QHFmFH = gbl_QHFmFH; %Q(~HF - F~H)
HWUsi= gbl_HWUsi;
W = gbl_W;

dUinv = {};
dU = {};
dUsqrtinv = {};
dY = {};
cIdY = {};
dHtilde = {};
dwGradE = {};

for k = [1:gbl_kpoints]
    dUinv{k} = -Uinv{k}*(dW{k}'*O(W{k}) + W{k}'*O(dW{k}))*Uinv{k};
    dU{k} = dW{k}'*O(W{k}) + W{k}'*O(dW{k});
    dUsqrtinv{k} = Q(dUinv{k}, Uinv{k});
    dY{k} = dW{k}*Usqrtinv{k}+W{k}*dUsqrtinv{k};
    cIdY{k} = cI(dY{k}, k);
end

dn = getdn(gbl_IY, cIdY, gbl_f);

excpn = gbl_excpn;
excppn = gbl_excppn;
JdagOJn = gbl_JdagOJn;
OJdN=O(cJ(dn));
dVsp = cJdag(O(-4*pi*Linv(OJdN))) ...
    + cJdag(O(cJ(excpn.*dn))) ...
    + Diagprod(excppn.*dn, JdagOJn) ...
    + Diagprod(excpn, cJdag(OJdN));

for k = [1:gbl_kpoints]
    dHW{k} = dH(IW{k}, dVsp, k);
end

for k = [1:gbl_kpoints]
    dHtilde{k} = (dUsqrtinv{k}*W{k}'+Usqrtinv{k}*dW{k}')*HWUsi{k};
    dHtilde{k} = dHtilde{k}+dHtilde{k}'+Usqrtinv{k}*(W{k}'*dHW{k})*Usqrtinv{k};
    
    tmp2{k} = dHW{k}*UsiFUsi{k} + H(dW{k}*UsiFUsi{k} + W{k}*dUsqrtinv{k}*F*Usqrtinv{k}+W{k}*Usqrtinv{k}*F*dUsqrtinv{k}, gbl_hsum, k);
    
    dwGradE{k} = -O(dW{k}*Uinv{k}*(W{k}'*tmp1{k}) + W{k}*dUinv{k}*(W{k}'*tmp1{k}) + W{k}*Uinv{k}*(dW{k}'*tmp1{k})) + tmp2{k} - O(W{k}*(Uinv{k}*(W{k}'*tmp2{k}))) ...
        + O(dY{k}*QHFmFH{k} + Y{k}*Q(dHtilde{k}*F-F*dHtilde{k}, U{k}) + Y{k}*dQ(HFmFH{k}, U{k}, dU{k}));
    dwGradE{k} = dwGradE{k}*gbl_weights(k);
end
end