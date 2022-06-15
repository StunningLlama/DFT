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
global gbl_IW;
global gbl_IY;
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

% T(k) k -> k, 1
% k T(k) -> Tk, 2
% T^-1(k) k -> k, 2
% k T^-1(k) -> Tk, 1

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dUinv{Tk_k} = -Uinv{k}*(dW{k_Tk}'*O(W{Tk}))*Uinv{Tk} - Uinv{Tk}*(W{Tk}'*O(dW{Tk_k}))*Uinv{k};
    dUinv{k_Tk} = -Uinv{T(k)}*(dW{Tk_k}'*O(W{k}))*Uinv{k} - Uinv{k}*(W{k}'*O(dW{Tk_k}))*Uinv{Tk};
end

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dU{Tk_k} = dW{k_Tk}'*O(W{Tk}) + W{Tk}'*O(dW{Tk_k});
    dU{k_Tk} = dW{Tk_k}'*O(W{k}) + W{k}'*O(dW{k_Tk});
end

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dUsqrtinv{Tk_k} = QInc(Uinv,dUinv,Tk_k);
    dUsqrtinv{k_Tk} = QInc(Uinv,dUinv,k_Tk);
end

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dY{Tk_k} = dW{Tk_k}*Usqrtinv{k}+W{Tk}*dUsqrtinv{Tk_k};
    dY{k_Tk} = dW{k_Tk}*Usqrtinv{Tk}+W{k}*dUsqrtinv{k_Tk};
end

dn = getdnInc(gbl_IY, cIdY, gbl_f);

excpn = gbl_excpn;
excppn = gbl_excppn;
JdagOJn = gbl_JdagOJn;
OJdN{q}=O(cJ(dn{q}));
OJdN{mq}=O(cJ(dn{mq}));
dVsp{q} = cJdag(O(-4*pi*Linv(OJdN{q}))) ...
    + cJdag(O(cJ(excpn.*dn{q}))) ...
    + Diagprod(excppn.*dn{q}, JdagOJn) ...
    + Diagprod(excpn, cJdag(OJdN{q}));
dVsp{mq} = cJdag(O(-4*pi*Linv(OJdN{mq}))) ...
    + cJdag(O(cJ(excpn.*dn{mq}))) ...
    + Diagprod(excppn.*dn{mq}, JdagOJn) ...
    + Diagprod(excpn, cJdag(OJdN{mq}));
%update q
k_Tk_q = gbl_kvectors(k)+gbl_kvectors(Tk)+gbl_kvectors(q);
Tk_mk_mq = gbl_kvectors(Tk)-gbl_kvectors(k)-gbl_kvectors(q);
for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dHW{Tk_k} = dHInc(gbl_IW{k}, M(dVsp{q}, k_Tk_q), k);
    dHW{k_Tk} = dHInc(gbl_IW{Tk}, M(dVsp{mq}, Tk_mk_mq), k);
end

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dHtilde{Tk_k} = dY{k_Tk}'*H(Y{Tk}, gbl_hsum, Tk) + Y{Tk}'*dHInc(gbl_IY{k}, M(dVsp{q}, k_Tk_q), k) + Y{Tk}'*H(dY{Tk_k}, gbl_hsum, Tk);
    dHtilde{k_Tk} = dY{Tk_k}'*H(Y{k}, gbl_hsum, k) + Y{k}'*dHInc(gbl_IY{Tk}, M(dVsp{q}, Tk_mk_mq), Tk) + Y{k}'*H(dY{k_Tk}, gbl_hsum, Tk);%todo
    tmp2{Tk_k} = dHW{Tk_k}*UsiFUsi{k} + H(dY{Tk_k}*F*Usqrtinv{k}+Y{Tk}*F*dUsqrtinv{Tk_k}, gbl_hsum, Tk);
    tmp2{k_Tk} = dHW{k_Tk}*UsiFUsi{Tk} + H(dY{k_Tk}*F*Usqrtinv{Tk}+Y{k}*F*dUsqrtinv{k_Tk}, gbl_hsum, Tk);
    
    dwGradE{Tk_k} = -O(dY{Tk_k}*(Y{k}'*tmp1{k})) - O(Y{k}*(dY{k_Tk}'*tmp1{Tk})) + tmp2{Tk_k} - O(W{Tk}*(Uinv{Tk}*(W{Tk}'*tmp2{Tk_k}))) ...
        + O(dY{Tk_k}*QHFmFH{k} + Y{Tk}*dQInc(Htilde, F, U, dU, Tk, k, Tk_k));
    dwGradE{Tk_k} = dwGradE{Tk_k}*gbl_weights(k);
    
    dwGradE{k_Tk} = -O(dY{k_Tk}*(Y{Tk}'*tmp1{Tk})) - O(Y{Tk}*(dY{Tk_k}'*tmp1{k})) + tmp2{k_Tk} - O(W{k}*(Uinv{k}*(W{k}'*tmp2{k_Tk}))) ...
        + O(dY{k_Tk}*QHFmFH{Tk} + Y{k}*dQInc(Htilde, F, U, dU, k, Tk, k_Tk));
    dwGradE{k_Tk} = dwGradE{k_Tk}*gbl_weights(k);
end
end