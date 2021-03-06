h=0.00000000001;
X=[8 8 8; 8+2 8 8];
kS_orig = [4; 2; 1];
comparen = true;

setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);
rng('default');
rng(246);
global gbl_kS;
q = coordtoindex([1 0 0], gbl_kS)+1;
W=initializeRandomState();
%W = sd(W, 5);
W = orthonormalize(W);
dWa=initializeRandomStateInc(q);

setupPccgWavefuncInc(W);
dGrad = getPsiPsiDerivWFillingsInc(dWa,q);

global TEST_A;
if (comparen)
    TEST_A{1} = cJcomp(TEST_A{1},q);
    TEST_A{2} = cJcomp(TEST_A{2},Tinv(1,q));
end

global gbl_f;
n = getn(W, gbl_f);
nn = {cJcomp(n,1)};

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);

Wbig = copyToLargeCell(W, kS_orig, q);
dWabig = copyToLargeCell(dWa, kS_orig, q);
dGradbig = copyToLargeCell(dGrad, kS_orig, q);
testbig = copyToLargeCell(TEST_A, kS_orig, q);
nbig = copyToLargeCell(nn, kS_orig, q);

Wbig{1} = Wbig{1}/sqrt(prod(kS_orig));
dWabig{1} = dWabig{1}/sqrt(prod(kS_orig));
dGradbig{1} = dGradbig{1}/sqrt(prod(kS_orig));
if (comparen)
    testbig = testbig;
end
setupPccgWavefunc(Wbig);
dGrad2 = getPsiPsiDerivWFillings(dWabig);
dGradFD2 = mult(linadd(getgrad(linadd(Wbig,dWabig,1,h)),getgrad(Wbig),1,-1),1/h);
n2 = getn(Wbig, gbl_f);

global TEST_B;

if (comparen)
    TEST_B = cJcomp(TEST_B,1);
end

if (comparen)
    RATIOtest = testbig./TEST_B;
else
    RATIOtest = testbig{1}./TEST_B{1};
end

RATIO =(dGrad2{1}./dGradbig{1});
nRATIO = nbig./cJcomp(n2,1);
%rat = dGradbig{1}./dGrad2{1};