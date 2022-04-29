function W = iterate(sdits)
W=initializeRandomState()

format long
%# Converge
W=sd(W,sdits); %# 20 iterations of simple sd() to get nearer to the minimum
W = orthonormalize(W); %# Restart as orthonormal functions
W=pccg(W,50,1); %# 50 iterations of pclm from same W
E = getE(W);
end