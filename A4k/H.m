function out=H(W, Vsp, k)
%Y=W*inv(sqrtm(W'*O(W)));
%n = getn(cI(Y), gbl_f);
%uinv = inv(W'*O(W));
%n = gbl_f*diagouter(cI(W*uinv), cI(W));
out = -0.5*L(W, k);
for col=1:size(W,2)
out(:,col) = out(:,col) + cIdag(Diagprod(Vsp, cI(W(:,col))));
end
end