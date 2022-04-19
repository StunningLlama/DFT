setupSmallGe();
ewald();
[W,E1] = iterate(20);
getE(W)
%global gbl_X;
%Visualize(W, X);
setupBigGe();
ewald();
[W,E1] = iterate(30);
getE(W)