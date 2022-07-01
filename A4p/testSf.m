
X=[8 8 8; 8+2 8 8];
setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), [1; 1; 1], false);

global gbl_G2; global gbl_X; global gbl_G; global gbl_Z;
G2 = gbl_G2;
X = gbl_X;
G = gbl_G;
Z = gbl_Z;

dX = randn(size(X))*0.0000001;


chargefactor = ones(size(G, 1), 1)*Z;
Sf=sum(chargefactor.*exp(-i*G*X'), 2);

Sf2=sum(chargefactor.*exp(-i*G*(X+dX)'), 2);

dSf = getdSf(X, dX);
disp2(Sf2-Sf)
disp2(dSf)