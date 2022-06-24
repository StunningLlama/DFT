function Wbig = copyToLargeCell(W, kS_orig, q)
global gbl_kS;
global gbl_activeindex;
global gbl_activelength;
global gbl_G2c;
kS_tmp = gbl_kS;
gbl_kS = kS_orig;
Wbig = [];
if (iscell(W))
    if (size(W,2)==2)
        mq = Tinv(1, q);
        Wbig = zeros(size(gbl_G2c{1}));
        Wbig([gbl_activeindex(q): gbl_activeindex(q)+gbl_activelength(q)-1],1) = W{1};
        Wbig([gbl_activeindex(mq): gbl_activeindex(mq)+gbl_activelength(mq)-1],1) = W{2};
    elseif (size(W,2) == prod(kS_orig))
        Wbig = initializeZeroState();
        for k = [1:prod(kS_orig)]
            Wbig{1} = setLargeCell(Wbig{1}, k, k, W{k});
        end
    elseif (size(W,2) == 2*prod(kS_orig))
        Wbig = initializeZeroState();
        for k = [1:prod(kS_orig)]
            Tk = T(k,q); Tk_k = k; k_Tk = k+prod(kS_orig);
            Wbig{1} = setLargeCell(Wbig{1}, Tk, k, W{Tk_k});
            Wbig{1} = setLargeCell(Wbig{1}, k, Tk, W{k_Tk});
        end
    end
end
gbl_kS = kS_tmp;
end