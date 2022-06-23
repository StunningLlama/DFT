%setupSmallGe(48);
X=[8 8 8; 8+2 8 8];
setup(X, 1, 1, [48; 48; 48], diag([16 16 16]), [4; 1; 1], false);
%W = iterate(20);
%W = orthonormalize(W);
setupPccgWavefunc(W);
global gbl_kS;
q = coordtoindex([1 0 0], gbl_kS);

dWa = initializeRandomState();
A = complexinnerprod(dWa, getPsiPsiDerivWFillings(dWa));