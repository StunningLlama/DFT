function gradout = getdGradtmp(W, dW)
dWin = dW(1:length(gbl_active),1) + i*dW(length(gbl_active):2*length(gbl_active),1);
grad = getPsiPsiDeriv(W, dWin);
gradout = [real(grad); imag(grad)];
end