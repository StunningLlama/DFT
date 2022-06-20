function kq = Tinv(k,q)
    global gbl_kS;
    kq = coordtoindex(mod(indextocoord(k, gbl_kS)-indextocoord(q, gbl_kS), gbl_kS), gbl_kS)+1;
end