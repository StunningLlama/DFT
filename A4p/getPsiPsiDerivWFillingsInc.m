% Calculates differential of grad E given differential dW.
function dwGradE = getPsiPsiDerivWFillingsInc(dW, q)
global gbl_f;
global gbl_Htilde;
global gbl_Y;
global gbl_tmp1;
global gbl_HFmFH;
global gbl_IY;
global gbl_kpoints;
global gbl_weights;
global gbl_kvectors;
global gbl_W;
global gbl_excpn;
global gbl_excppn;
global gbl_JdagOJn;
global gbl_hsum;
global gbl_output;
F = diag(gbl_f);

Htilde = gbl_Htilde;
Y = gbl_Y;
tmp1 = gbl_tmp1;
HFmFH = gbl_HFmFH; %~HF - F~H
W = gbl_W;

dUinv = {};
dU = {};
dUsqrtinv = {};
dY = {};
cIdY = {};
dHtilde = {};
dwGradE = {};

% T(k) k -> k, 1
% k T(k) -> Tk, 2
% T^-1(k) k -> k, 2
% k T^-1(k) -> Tk, 1

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dU{Tk_k} = dW{k_Tk}'*O(W{k}) + W{Tk}'*O(dW{Tk_k});
    dU{k_Tk} = dW{Tk_k}'*O(W{Tk}) + W{k}'*O(dW{k_Tk});
    
    dUinv{Tk_k} = -dW{k_Tk}'*O(W{k}) - W{Tk}'*O(dW{Tk_k});
    dUinv{k_Tk} = -dW{Tk_k}'*O(W{Tk}) - W{k}'*O(dW{k_Tk});
    
    dUsqrtinv{Tk_k} = dUinv{Tk_k}/2.0;
    dUsqrtinv{k_Tk} = dUinv{k_Tk}/2.0;
end

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dY{Tk_k} = dW{Tk_k}+W{Tk}*dUsqrtinv{Tk_k};
    dY{k_Tk} = dW{k_Tk}+W{k}*dUsqrtinv{k_Tk};
    
    cIdY{Tk_k} = cI(dY{Tk_k}, Tk);
    cIdY{k_Tk} = cI(dY{k_Tk}, k);
end

global gbl_dY;
gbl_dY = dY;

dn = getdnInc(gbl_IY, cIdY, gbl_f, q);
gbl_output = dn;

qc = 1;
mqc = 2;
mq = Tinv(1, q);

excpn = gbl_excpn;
excppn = gbl_excppn;
JdagOJn = gbl_JdagOJn;
OJdN{qc}=O(cJ(dn{qc}));
OJdN{mqc}=O(cJ(dn{mqc}));
dVsp{qc} = cJdag(O(-4*pi*LinvInc(OJdN{qc},q,1))) ...
    + cJdag(O(cJ(excpn.*dn{qc}))) ...
    + Diagprod(excppn.*dn{qc}, JdagOJn) ...
    + Diagprod(excpn, cJdag(OJdN{qc}));
dVsp{mqc} = cJdag(O(-4*pi*LinvInc(OJdN{mqc},mq,1))) ...
    + cJdag(O(cJ(excpn.*dn{mqc}))) ...
    + Diagprod(excppn.*dn{mqc}, JdagOJn) ...
    + Diagprod(excpn, cJdag(OJdN{mqc}));

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    k_mTk_q = gbl_kvectors(k,:)-gbl_kvectors(Tk,:)+gbl_kvectors(q,:);
    Tk_mk_mq = gbl_kvectors(Tk,:)-gbl_kvectors(k,:)+gbl_kvectors(mq,:);
    
    dHY{Tk_k} = dH(gbl_IY{k}, M(k_mTk_q, dVsp{qc}), Tk);
    dHY{k_Tk} = dH(gbl_IY{Tk}, M(Tk_mk_mq, dVsp{mqc}), k);
end

global TEST_A;
TEST_A = dVsp;

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    
    dHtilde{Tk_k} = dY{k_Tk}'*H(Y{k}, gbl_hsum, k) + Y{Tk}'*dHY{Tk_k} + Y{Tk}'*H(dY{Tk_k}, gbl_hsum, Tk);
    dHtilde{k_Tk} = dY{Tk_k}'*H(Y{Tk}, gbl_hsum, Tk) + Y{k}'*dHY{k_Tk} + Y{k}'*H(dY{k_Tk}, gbl_hsum, k);
    tmp2{Tk_k} = dHY{Tk_k}*F + H(dY{Tk_k}*F+Y{Tk}*F*dUsqrtinv{Tk_k}, gbl_hsum, Tk);
    tmp2{k_Tk} = dHY{k_Tk}*F + H(dY{k_Tk}*F+Y{k}*F*dUsqrtinv{k_Tk}, gbl_hsum, k);
end

%TODO test positive definite operator?

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dwGradE{Tk_k} = -O(dY{Tk_k}*(Y{k}'*tmp1{k}) + Y{Tk}*(dY{k_Tk}'*tmp1{k})) + tmp2{Tk_k} - O(W{Tk}*(W{Tk}'*tmp2{Tk_k})) ...
        + O(dY{Tk_k}*HFmFH{k}/2.0 + Y{Tk}*((dHtilde{Tk_k}*F-F*dHtilde{Tk_k})/2.0 - (HFmFH{Tk}*dU{Tk_k} + dU{Tk_k}*HFmFH{k})/8.0));
    dwGradE{Tk_k} = dwGradE{Tk_k}*gbl_weights(k);
    
    dwGradE{k_Tk} = -O(dY{k_Tk}*(Y{Tk}'*tmp1{Tk}) + Y{k}*(dY{Tk_k}'*tmp1{Tk})) + tmp2{k_Tk} - O(W{k}*(W{k}'*tmp2{k_Tk})) ...
        + O(dY{k_Tk}*HFmFH{Tk}/2.0 + Y{k}*((dHtilde{k_Tk}*F-F*dHtilde{k_Tk})/2.0 - (HFmFH{k}*dU{k_Tk} + dU{k_Tk}*HFmFH{Tk})/8.0));
    dwGradE{k_Tk} = dwGradE{k_Tk}*gbl_weights(k);
end
end