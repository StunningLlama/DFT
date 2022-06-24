function kq = Tinv(k,q)
    global gbl_kS;
    kq = coordtoindex(mod(indextocoord(k-1, gbl_kS)-indextocoord(q-1, gbl_kS), gbl_kS'), gbl_kS)+1;
end