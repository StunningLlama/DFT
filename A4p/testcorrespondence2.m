h=0.00000000001;
X=[8 8 8; 8+2 8 8];
kS_orig = [4; 1; 1];
comparen = true;

setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);
rng('default');
rng(246);
global gbl_kS;
q = coordtoindex([1 0 0], gbl_kS)+1;
W=initializeRandomState();
W = sd(W, 5);
W = orthonormalize(W);
dWa=initializeRandomStateInc(q);

setupPccgWavefuncInc(W);
dGrad = getPsiPsiDerivWFillingsInc(dWa,q);

global TEST_A;
if (comparen)
    TEST_A{1} = cIdag(TEST_A{1},q);
    TEST_A{2} = cIdag(TEST_A{2},Tinv(1,q));
end

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);
Wbig = initializeZeroState();
dWabig = initializeZeroState();
dGradbig = initializeZeroState();
testbig = initializeZeroState();

Wbig = copyToLargeCell(W, kS_orig, q);
dWabig = copyToLargeCell(dWa, kS_orig, q);
dGradbig = copyToLargeCell(dGrad, kS_orig, q);
testbig = copyToLargeCell(TEST_A, kS_orig, q);

Wbig{1} = Wbig{1}/sqrt(prod(kS_orig));
dWabig{1} = dWabig{1}/(prod(kS_orig));
dGradbig{1} = dGradbig{1}*prod(kS_orig);
if (comparen)
    testbig = testbig*sqrt(prod(kS_orig));
end
setupPccgWavefunc(Wbig);
dGrad2 = getPsiPsiDerivWFillings(dWabig);
dGradFD2 = mult(linadd(getgrad(linadd(Wbig,dWabig,1,h)),getgrad(Wbig),1,-1),1/h);

global TEST_B;

if (comparen)
    TEST_B = cIdag(TEST_B,1);
end

if (comparen)
    RATIOtest = testbig./TEST_B;
else
    RATIOtest = testbig{1}./TEST_B{1};
end

RATIO =(dGrad2{1}./dGradbig{1});
%rat = dGradbig{1}./dGrad2{1};