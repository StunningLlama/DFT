function k = calcSpringConstant(W, X, dXa, dXb, dWa, dWb)
E1 = getTauTauDeriv(W, X, dXa, dXb);
%E2 = 2*real(trace(dWa'*getPsiTauDerivWFillings(W, dXb)));
E3 = 2*real(trace(dWb'*getPsiTauDerivWFillings(W, dXa)))
%E4 = 2*real(trace(dWa'*getPsiPsiDerivWFillings(W, dWb)))
%k = -(E1+E2+E3+E4);
k = -(E1+E3);
end