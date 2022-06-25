% Calculates differential of grad E given differential dW.
function dwGradE = getPsiPsiDerivWFillings(dW)
global gbl_f;
global gbl_Y;
global gbl_tmp1;
global gbl_IY;
global gbl_kpoints;
global gbl_weights;
global gbl_W;
global gbl_excpn;
global gbl_excppn;
global gbl_JdagOJn;
global gbl_hsum;

global gbl_dY;
global gbl_output;

F = diag(gbl_f);

Y = gbl_Y;
IY = gbl_IY;
tmp1 = gbl_tmp1;
W = gbl_W;

dUinv = {};
dU = {};
dUsqrtinv = {};
dY = {};
cIdY = {};
dHtilde = {};
dwGradE = {};

for k = [1:gbl_kpoints]
    dU{k} = dW{k}'*O(W{k}) + W{k}'*O(dW{k});
    dUinv{k} = -dU{k};
    dUsqrtinv{k} = dUinv{k}/2.0;
    dY{k} = dW{k}+W{k}*dUsqrtinv{k};
    cIdY{k} = cI(dY{k}, k);
end

gbl_dY = dY;

dn = getdn(gbl_IY, cIdY, gbl_f);
gbl_output = dn;

excpn = gbl_excpn;
excppn = gbl_excppn;
JdagOJn = gbl_JdagOJn;
OJdN=O(cJ(dn));
dVsp = cJdag(O(-4*pi*Linv(OJdN))) ...
    + cJdag(O(cJ(excpn.*dn))) ...
    + Diagprod(excppn.*dn, JdagOJn) ...
    + Diagprod(excpn, cJdag(OJdN));

for k = [1:gbl_kpoints]
    dHY{k} = dH(IY{k}, dVsp, k);
end
global TEST_B;
TEST_B = dVsp;

for k = [1:gbl_kpoints]
    dHtilde{k} = dY{k}'*H(Y{k},gbl_hsum,k) + Y{k}'*dH(IY{k}, dVsp, k) + Y{k}'*H(dY{k},gbl_hsum,k);
    
    tmp2{k} = dH(IY{k}, dVsp, k)*F + H(dY{k}*F + Y{k}*F*dUsqrtinv{k},gbl_hsum,k);
    
    dwGradE{k} = -O(dY{k}*(Y{k}'*tmp1{k}) + Y{k}*(dY{k}'*tmp1{k})) + tmp2{k} - O(Y{k}*(Y{k}'*tmp2{k})); %...
        %+ O(dY{k}*QHFmFH{k} + Y{k}*Q(dHtilde{k}*F-F*dHtilde{k}, U{k}) + Y{k}*dQ(HFmFH{k}, U{k}, dU{k}));
    dwGradE{k} = dwGradE{k}*gbl_weights(k);
end
end