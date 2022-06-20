function dn = getdnInc(IY, IdY, f, q)
global gbl_G2;
dn = {};
qc = 1;
mqc = 2;
dn{qc} = zeros(size(gbl_G2));
dn{mqc} = zeros(size(gbl_G2));
global gbl_kpoints;
global gbl_weights;
global gbl_kvectors;
for k = [1:gbl_kpoints]
    Tk = T(k,q); Tinvk = Tinv(k,q); Tk_k = k; Tinvk_k = Tinv(k,q)+gbl_kpoints;
    %Tinv(k) k = Tinv(k), T Tinv(k)
    
    fillings = (ones(size(IY{k}, 1), 1)*f');
    
    k_mTinvk_mq = gbl_kvectors(k)+gbl_kvectors(Tinvk)+gbl_kvectors(qc);
    mk_Tk_mq = gbl_kvectors(k)+gbl_kvectors(Tk)+gbl_kvectors(qc);
    k_mTk_q = gbl_kvectors(k)+gbl_kvectors(Tk)+gbl_kvectors(qc);
    mk_Tinvk_q = gbl_kvectors(k)+gbl_kvectors(Tinvk)+gbl_kvectors(qc);
    dn{qc} = dn{qc} + M(k_mTinvk_mq, sum(fillings.*conj(IdY{Tinvk_k}).*IY{k},2))+ ...
        M(mk_Tk_mq, sum(fillings.*conj(IY{k}).*IdY{Tk_k}, 2));
    dn{mqc} = dn{mqc} + M(k_mTk_q, sum(fillings.*conj(IdY{Tk_k}).*IY{k},2))+ ...
        M(mk_Tinvk_q, sum(fillings.*conj(IY{k}).*IdY{Tinvk_k}, 2));
end
end