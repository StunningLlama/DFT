setupSmallGe(16);
ewald()
W = iterate(20);
E0 = getE(W)
disp("Eigenvalues");
[Psi, epsilon] = getPsi(W);
epsilon

%global gbl_X;
%Visualize(W, X);

setupBigGe(16);
ewald()
W = iterate(30);
E1 = getE(W)
disp("Eigenvalues");
[Psi, epsilon] = getPsi(W);
epsilon
E1/E0