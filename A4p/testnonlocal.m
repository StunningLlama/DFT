
X=[8 8 8; 8+2 8 8];
setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), [1;1;1], false, true);
global gbl_Vnl;
%visualize2(real(gbl_Vnl));
W = iterate(20);
visualize(W);
