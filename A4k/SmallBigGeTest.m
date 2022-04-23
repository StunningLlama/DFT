setupSmallGe(48);
ewald()
[W,E] = iterate(20);
getE(W)
disp("Eigenvalues");
[Psi, epsilon] = getPsi(W);
epsilon

%global gbl_X;
%Visualize(W, X);

setupBigGe(48);
ewald()
[W,E] = iterate(30);
getE(W)
disp("Eigenvalues");
[Psi, epislon] = getPsi(W);
epsilon