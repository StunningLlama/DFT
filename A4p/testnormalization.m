global gbl_R; global gbl_kS; global gbl_r; global gbl_output; global gbl_kvectors;

X=[8 8 8; 8+2 8 8];
kS_orig = [4; 1; 1];

setup(X, 1, [1 1], [32; 32; 32], diag([16 16 16]), kS_orig, false);
global gbl_G2c;
test = zeros(size(gbl_G2c{1}));
test(1)=1;
A = cI(test,1);

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);
testt = zeros(size(gbl_G2c{1}));
testt(1)=1;
B = cI(testt,1);