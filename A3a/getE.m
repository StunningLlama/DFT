function E=getE(W)
global gbl_Vdual;
global gbl_f;
uinv = inv(W'*O(W));
%n = gbl_f*diagouter(cI(W*uinv), cI(W));
Y=W*sqrtm(uinv);
n = getn(cI(Y), gbl_f);
E = real(-0.5*gbl_f*trace(W'*L(W)*uinv) + gbl_Vdual'*n + 0.5*n'*cJdag(O(-4*pi*Linv(O(cJ(n))))) + n'*cJdag(O(cJ(excVWN(n)))));
end