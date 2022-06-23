X=[8 8 8; 8+2 8 8];
% setup(X, 1, 1, [48; 48; 48], diag([16 16 16]), [2; 2; 2], false);
% W = iterate(20);
% getE(W)

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), [2; 2; 2], false);
W = iterate(20);
getE(W)