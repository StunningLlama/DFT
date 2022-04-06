
global gbl_G2; global gbl_X; global gbl_G; global gbl_Z;
G2 = gbl_G2;
X = gbl_X;
G = gbl_G;
Z = gbl_Z;

X = randn(4, 3);
dX = randn(4, 3)*0.0000001;


chargefactor = ones(size(G, 1), 1)*Z;
Sf=sum(chargefactor.*exp(-i*G*X'), 2);

Sf2=sum(chargefactor.*exp(-i*G*(X+dX)'), 2);

dSf = getdSf(X, dX);
disp2(Sf2-Sf)
disp2(dSf)