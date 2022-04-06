function out=H2(W, dW)
global gbl_hsum;


%uinv = inv(W'*O(W));
%n = gbl_f*diagouter(cI(W*uinv), cI(W));
out = -0.5*L(dW);
for col=1:size(W,2)
out(:,col) = out(:,col) + cIdag(Diagprod(gbl_hsum, cI(dW(:,col))));
end
end