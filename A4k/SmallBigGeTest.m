setupSmallGe();
ewald();
[W,E1] = iterate(20);
getE(W)
disp("Eigenvalues");
[Psi, epsilon] = getPsi(W);
epsilon
%global gbl_X;
%Visualize(W, X);
setupBigGe();
ewald();
[W,E1] = iterate(30);
getE(W)
disp("Eigenvalues");
[Psi, epislon] = getPsi(W);
epsilon