% Performs finite difference test of getE() and getgrad()
%
% Usage: fdtest(W,S)
%
% W: starting point for test (size: prod(S) x Ns)
% S: Dimensions of 3d data
function fdtest()
%setupSmallGe(32);
setup([8 8 8; 8+2 8 8], 1, [1 1], [48; 48; 48], diag([16 16 16]), [1;1;1], false, true);

W = initializeRandomState();
%# Compute intial energy and gradient
E0=getE(W)
g0=getgrad(W);
%# Choose a random direction to explore
dW=initializeRandomState();
%# Explore a range of step sizes decreasing by powers of ten
for delta=10.^[1:-1:-9]
%# Directional derivative formula
dE=delta*2*complexinnerprod(g0, dW);
%# Print ratio of actual change to expected change, along with estimate
%# of the error in this quantity due to rounding
fprintf(' %20.16f\n %20.16f\n\n', ...
(getE(linadd(W,dW,1,delta))-E0)/dE, sqrt(size(W{1},1))*eps/abs(dE) );
end
end