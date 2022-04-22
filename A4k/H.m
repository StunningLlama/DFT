function out=H(W, n, k)
global gbl_Vdual;
global gbl_f;
%Y=W*inv(sqrtm(W'*O(W)));
%n = getn(cI(Y), gbl_f);
%uinv = inv(W'*O(W));
%n = gbl_f*diagouter(cI(W*uinv), cI(W));
out = -0.5*L(W, k);
lhs = real(gbl_Vdual) + cJdag(O(-4*pi*Linv(O(cJ(n))))) ...
    + cJdag(O(cJ(excVWN(n)))) ...
    + Diagprod(excpVWN(n), cJdag(O(cJ(n))));
for col=1:size(W,2)
out(:,col) = out(:,col) + cIdag(Diagprod(lhs, cI(W(:,col))));
end
end