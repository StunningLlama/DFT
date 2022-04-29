function k = calcSpringConstant(W, X, dXa, dXb, dWa, dWb)
E1 = getTauTauDeriv(dXa, dXb);
%E2 = 2*real(trace(dWa'*getPsiTauDerivWFillings(dXb)));
E3 = 2*complexinnerprod(dWb, getPsiTauDerivWFillings(dXa));
%E3b = 2*real(trace(dWb'*getPsiTauDerivWFillings(dXa)))
%E4 = 2*real(trace(dWa'*getPsiPsiDerivWFillings(dWb)));
k = -(E1+E3);
end