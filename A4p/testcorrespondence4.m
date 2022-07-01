h=0.00000001;
global gbl_Vtest;
gbl_Vtest = 0;
X=[8 8 8; 8+2 8 8];
kS_orig = [4; 2; 1];

setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);
rng('default');
rng(246);
global gbl_kS; global gbl_R; global gbl_r; global gbl_kvectors;
q = coordtoindex([1 0 0], gbl_kS)+1;

pk = [1.0 0 0]*2*pi*inv(gbl_R);
p = sum(gbl_r*diag(pk),2);
dVsp = {};
dVsp{1} = sin(p);
dVsp{2} = M(-gbl_kvectors(q,:)-gbl_kvectors(Tinv(1,q),:), sin(p));
cJdVsp = {};
cJdVsp{1} = cJcomp(dVsp{1},q);
cJdVsp{2} = cJcomp(dVsp{2},Tinv(1,q));

W=initializeRandomState();
%W = sd(W, 5);
W = orthonormalize(W);

setupPccgWavefuncInc(W);
dGrad = getPsiVspDerivWFillingsInc(dVsp,q);

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);


dVspBig = cI(copyToLargeCell(cJdVsp, kS_orig, q), 1);

Wbig = copyToLargeCell(W, kS_orig, q);
dGradbig = copyToLargeCell(dGrad, kS_orig, q);

Wbig{1} = Wbig{1}/sqrt(prod(kS_orig));
dGradbig{1} = dGradbig{1};

setupPccgWavefunc(Wbig);
dGrad2 = getPsiVspDerivWFillings(dVspBig);
G1 = getgrad(Wbig);
global gbl_Vdual;
gbl_Vdual = gbl_Vdual + dVspBig*h;
G2 = getgrad(Wbig);
dGradFD2 = mult(linadd(G2,G1,1,-1),1/h);

RATIOFD =(dGrad2{1}./dGradFD2{1});
RATIO =(dGrad2{1}./dGradbig{1});