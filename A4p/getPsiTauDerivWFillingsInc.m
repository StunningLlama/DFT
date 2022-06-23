% Calculates differential of grad E given differential dTau.
function dwGradE = getPsiTauDerivWFillingsInc(dVsp, q)

global gbl_f;
global gbl_Y;
global gbl_IW;
global gbl_tmp1;
global gbl_kpoints;
global gbl_weights;
global gbl_kvectors;

qc = 1;
mqc = 2;

F = diag(gbl_f);

Y = gbl_Y;
tmp1 = gbl_tmp1;

dHtilde = {};
dwGradE = {};

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    k_Tk_q = gbl_kvectors(k,:)+gbl_kvectors(Tk,:)+gbl_kvectors(q,:);
    Tk_mk_mq = gbl_kvectors(Tk,:)-gbl_kvectors(k,:)-gbl_kvectors(q,:);
    
    dHY{Tk_k} = dH(gbl_IW{k}, M(k_Tk_q, dVsp{qc}), Tk);
    dHY{k_Tk} = dH(gbl_IW{Tk}, M(Tk_mk_mq, dVsp{mqc}), k);
end

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    dHtilde{Tk_k} = Y{Tk}'*dHY{Tk_k};
    dHtilde{k_Tk} = Y{k}'*dHY{k_Tk};
end

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tk_k = k; k_Tk = k+gbl_kpoints;
    tmp1{Tk_k} = dHY{Tk_k}*F;
    tmp1{k_Tk} = dHY{k_Tk}*F;
    
    dwGradE{Tk_k} = tmp1{Tk_k} - O(Y{Tk}*Y{Tk}'*tmp1{Tk_k}) + O(Y{Tk}*((dHtilde{Tk_k}*F-F*dHtilde{Tk_k})/2.0));
    dwGradE{Tk_k} = dwGradE{Tk_k}*gbl_weights(k);
    
    dwGradE{k_Tk} = tmp1{k_Tk} - O(Y{k}*Y{k}'*tmp1{k_Tk}) + O(Y{k}*((dHtilde{k_Tk}*F-F*dHtilde{k_Tk})/2.0));
    dwGradE{k_Tk} = dwGradE{k_Tk}*gbl_weights(k);
end
end