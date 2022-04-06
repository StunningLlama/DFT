function Sf = getSf(X)
global gbl_Z;
global gbl_G;
Z = gbl_Z;
G = gbl_G;
chargefactor = ones(size(G, 1), 1)*Z;
Sf=sum(chargefactor.*exp(-i*G*X'), 2);
end