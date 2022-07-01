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
dVdual=cJdag(dVtilde);


global gbl_K;
global gbl_Mnl;
global gbl_Vnl;
global gbl_Z;
global gbl_Gc;

if (k==1)
    chargefactor = ones(size(gbl_Gc{1}, 1), 1)*gbl_Z;
    dSf2=sum(-i*(gbl_Gc{1}*dX').*chargefactor.*exp(-i*gbl_Gc{1}*X'), 2);
    dK = O(cJcomp(gbl_Vnl, 1).*dSf2);
    
    out = dK*(gbl_Mnl*(gbl_K'*W)) + gbl_K*(gbl_Mnl*(dK'*W));
else
    out = zeros(size(W));
end

for col=1:size(W,2)
out(:,col) = out(:,col) + cIdag(Diagprod(dVdual, cI(W(:,col), k)), k);
end
end