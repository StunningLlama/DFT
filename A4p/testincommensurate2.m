global gbl_R; global gbl_kS; global gbl_r; global gbl_output; global gbl_kvectors;

X=[8 8 8; 8+2 8 8];
kS_orig = [4; 2; 1];

setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);
W = iterate(20);
W = orthonormalize(W);

q = coordtoindex([1 0 0], gbl_kS)+1;
pk = [1.0 0 0]*2*pi*inv(gbl_R);
p = sum(gbl_r*diag(pk),2);
dVsp = {};
dVsp{1} = sin(p);
dVsp{2} = M(-gbl_kvectors(q,:)-gbl_kvectors(Tinv(1,q),:), sin(p));
%visualize2(dVsp{1});
%visualize2(-dVsp{1});
setupPccgWavefuncInc(W);
pccgWavefuncInc(W, {cJdag(O(cJ(dVsp{1}))),cJdag(O(cJ(dVsp{2})))}, q, 100, 1);
dn = gbl_output;

cJdVsp = {};
cJdVsp{1} = cJcomp(dVsp{1},q);
cJdVsp{2} = cJcomp(dVsp{2},Tinv(1,q));
cJdn = {};
cJdn{1} = cJcomp(dn{1},q);
cJdn{2} = cJcomp(dn{2},Tinv(1,q));

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), kS_orig, false);

Wbig = copyToLargeCell(W, kS_orig, q);
dVspBig = cI(copyToLargeCell(cJdVsp, kS_orig, q), 1);
dnnew = copyToLargeCell(cJdn, kS_orig, q);

Wbig{1} = Wbig{1}/sqrt(prod(kS_orig));

%visualize2(real(dVspBig));
%visualize2(-real(dVspBig));
setupPccgWavefunc(Wbig);
pccgWavefuncVsp(Wbig, cJdag(O(cJ(dVspBig))), 100, 1);
dnold = cJcomp(gbl_output, 1);
A = dnnew./dnold;