function d2U = getdsqEwald(X, dXa, dXb)
global gbl_S; global gbl_R; global gbl_G2; global gbl_X; global gbl_G; global gbl_Sf; global gbl_X; global gbl_f; global gbl_r; global gbl_Z;

S = gbl_S;
R = gbl_R;
G2 = gbl_G2;
X = gbl_X;
G = gbl_G;
Sf = gbl_Sf;
X = gbl_X;
r = gbl_r;
Z = gbl_Z;

dr= sqrt(sum((r - ones(prod(S), 1)*sum(R,2)'/2).^2, 2));
sigma1=0.25;
g1=exp(-dr.^2/(2*sigma1^2))/sqrt(2*pi*sigma1^2)^3;

daSf = getdSf(X, dXa);
dbSf = getdSf(X, dXb);
d2Sf = getdsqSf(X, dXa, dXb);

nnuc=real(cI(cJ(g1).*Sf));
phi=real(cI(Linv(-4*pi*O(cJ(nnuc)))));


dannuc=real(cI(cJ(g1).*daSf));
dbnnuc=real(cI(cJ(g1).*dbSf));

daphi=real(cI(Linv(-4*pi*O(cJ(dannuc)))));
dbphi=real(cI(Linv(-4*pi*O(cJ(dbnnuc)))));

d2nnuc=real(cI(cJ(g1).*d2Sf));
d2phi=real(cI(Linv(-4*pi*O(cJ(d2nnuc)))));
d2U = 0.5*real(cJ(d2phi)'*O(cJ(nnuc)) +...
    cJ(daphi)'*O(cJ(dbnnuc)) +...
    cJ(dbphi)'*O(cJ(dannuc)) +...
    cJ(phi)'*O(cJ(d2nnuc)));
end