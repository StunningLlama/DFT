function out=KInc(W, q)
global gbl_Gc;
global gbl_R;
global gbl_kvectors;
global gbl_weights;
global gbl_kpoints;
out = {};

for k = [1:gbl_kpoints]
    Tk = T(k,q); Tinvk = Tinv(k,q); Tk_k = k; Tinvk_k = Tinv(k,q)+gbl_kpoints;
    
    Tkvec = gbl_kvectors(Tk,:);
    Tinvkvec = gbl_kvectors(Tinvk,:);
    karray_Tk_k = ones(size(W{Tk_k}, 1), 1)*Tkvec;
    karray_Tinvk_k = ones(size(W{Tinvk_k}, 1), 1)*Tinvkvec;
    
    %out{k} = -(1/gbl_weights(k))*W{k}./(det(gbl_R)*((sum((gbl_Gc+karray).^2, 2)+1)*ones(1,size(W,2))));
    out{Tk_k} = (1/gbl_weights(k))*W{Tk_k}./((sum((gbl_Gc{Tk}+karray_Tk_k).^2, 2)+1)*ones(1,size(W{Tk_k},2)));
    out{Tinvk_k} = (1/gbl_weights(k))*W{Tinvk_k}./((sum((gbl_Gc{Tinvk}+karray_Tinvk_k).^2, 2)+1)*ones(1,size(W{Tinvk_k},2)));
end
end