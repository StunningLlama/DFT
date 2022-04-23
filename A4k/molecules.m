molecule = 2s
if (molecule==1) % Dihelium
    X=[6 6 9; 6+1.5 6 9];
    setup(X, 1, [1 1], true);
    [W,E1] = iterate(20);
    visualize(W, X);
elseif (molecule==2) % Dihelium
    X=[8 8 8; 8+1.5 8 8];
    setup(X, 2, [2 2], true);
    [W,E1] = iterate(20);
    visualize(W, X);
elseif (molecule==3) % Methane
    bondlen = 108.7*0.0188972599;
    angle0 = deg2rad(109.5);
    angle1 = deg2rad(0);
    angle2 = deg2rad(120);
    angle3 = deg2rad(240);
    X=[8 8 8; 8 8 8; 8 8 8; 8 8 8; 8 8 8] + bondlen*[0 0 0; cos(angle1) sin(angle1) cos(angle0); cos(angle2) sin(angle2) cos(angle0); cos(angle3) sin(angle3) cos(angle0); 0 0 1];
    Z = [6 1 1 1 1];
    nstates = 5;
    setup(X, nstates, Z, true);
    [W,E1] = iterate(20);
    visualize(W, X);
elseif (molecule==4) % Crystal
    X=[1.5 1.5 1.5; 1.5+1 1.5 1.5];
    setup(X, 1, [1 1], false);
    [W,E1] = iterate(30);
    visualize(W, X);
end