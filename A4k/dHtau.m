function out = dHtau(W, dX, k)
global gbl_Vps;

global gbl_X;
X = gbl_X;

%Y=W*inv(sqrtm(W'*O(W)));
%n = getn(cI(Y), gbl_f);
%uinv = inv(W'*O(W));
%n = gbl_f*diagouter(cI(W*uinv), cI(W));


dSf = getdSf(X, dX);
dVtilde = gbl_Vps.*dSf;
dVtilde(1)=0.;
dVdual=cJ(dVtilde);


out = zeros(size(W));
for col=1:size(W,2)
out(:,col) = out(:,col) + cIdag(Diagprod(dVdual, cI(W(:,col), k)), k);
end
end