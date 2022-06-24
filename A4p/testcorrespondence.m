h=0.00000000001;
X=[8 8 8; 8+2 8 8];
setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), [2; 2; 2], false);
rng('default');
rng(246);
W=initializeRandomState();
W = sd(W, 5);
W = orthonormalize(W);
dWa=initializeRandomState();

setupPccgWavefunc(W);
dGrad = getPsiPsiDerivWFillings(dWa);
dGradFD = mult(linadd(getgrad(linadd(W,dWa,1,h)),getgrad(W),1,-1),1/h);

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), [2; 2; 2], false);
Wbig = initializeZeroState();
dWabig = initializeZeroState();
dGradbig = initializeZeroState();
dGradFDbig = initializeZeroState();
for k = [1:8]
    Wbig{1} = setLargeCell(Wbig{1}, k, k, W{k});
    dWabig{1} = setLargeCell(dWabig{1}, k, k, dWa{k});
    dGradbig{1} = setLargeCell(dGradbig{1}, k, k, dGrad{k})
    dGradFDbig{1} = setLargeCell(dGradFDbig{1}, k, k, dGradFD{k})
end
Wbig{1} = Wbig{1}/sqrt(8);
dWabig{1} = dWabig{1}/8;
setupPccgWavefunc(Wbig);
dGrad2 = getPsiPsiDerivWFillings(dWabig);
dGradFD2 = mult(linadd(getgrad(linadd(Wbig,dWabig,1,h)),getgrad(Wbig),1,-1),1/h);
A=(dGrad2{1}./dGradbig{1});
B=(dGradFD2{1}./dGradbig{1});
%rat = dGradbig{1}./dGrad2{1};