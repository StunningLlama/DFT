function dSf = getdSf(X, dX)
global gbl_Z;
global gbl_G;
Z = gbl_Z;
G = gbl_G;
chargefactor = ones(size(G, 1), 1)*Z;
dSf = sum(-i*(G*dX').*chargefactor.*(exp(-i*G*X')), 2);
end