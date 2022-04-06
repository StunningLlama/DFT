molecule = 1
if (molecule==1) % Dihelium
    X=[8 8 8; 8+1.5 8 8];
    setup(X, 1, [1 1]);
    [W,E1] = iterate();
    visualize(W);
elseif (molecule==2) % Dihelium
    X=[8 8 8; 8+1.5 8 8];
    setup(X, 2, [2 2]);
    [W,E1] = iterate();
    visualize(W);
elseif (molecule==3) % Methane
    X=[8 8 8; 8 8 8; 8 8 8; 8 8 8; 8 8 8] + [0 0 0; 0 1.5 -0.3; 1 -0.5 -0.3; -1 -0.5 -0.3; 0 1.5 0];
    setup(X, 5, [6 1 1 1 1]);
    [W,E1] = iterate();
    visualize(W);
end