function E=getE(W)
global gbl_Vdual;
uinv = inv(W'*O(W));
E = real(-0.5*trace(W'*L(W)*uinv) + gbl_Vdual'*diagouter(cI(W*uinv), cI(W)));
end