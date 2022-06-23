%setupSmallGe(48);
X=[8 8 8; 8+2 8 8];
setup(X, 1, 1, [48; 48; 48], diag([16 16 16]), [4; 1; 1], false);
W = iterate(20);
W = orthonormalize(W);
setupPccgWavefuncInc(W);
q = coordtoindex([1 0 0], gbl_kS);

dwA = initializeRandomStateInc(q);
dwB = initializeRandomStateInc(q);

A = getPsiPsiDerivWFillingsInc(linadd(dwA,dwB,1,1), q);
B = linadd(getPsiPsiDerivWFillingsInc(dwA, q), getPsiPsiDerivWFillingsInc(dwB, q), 1, 1);