function setupSmallGe(resolution)

S=[resolution; resolution; resolution];
a=5.66/0.52917721; %# Lattice constant (converted from angstroms to bohrs)
R=a*diag(ones(3,1));

X=a*[0.00 0.00 0.00 %# diamond lattice in cubic cell
0.25 0.25 0.25
0.00 0.50 0.50
0.25 0.75 0.75
0.50 0.00 0.50
0.75 0.25 0.75
0.50 0.50 0.00
0.75 0.75 0.25];
Z = 4*ones(1, 8);

Ns=16;

setup(X, Ns, Z, S, R, [2; 2; 2], true);
end