function setupBigGe(resolution)

S=[resolution*2; resolution*2; resolution*2];
a=5.66/0.52917721; %# Lattice constant (converted from angstroms to bohrs)
R=2*a*diag(ones(3,1));

X=[0.00 0.00 0.00 %# diamond lattice in cubic cell
0.25 0.25 0.25
0.00 0.50 0.50
0.25 0.75 0.75
0.50 0.00 0.50
0.75 0.25 0.75
0.50 0.50 0.00
0.75 0.75 0.25];

X = [X; X+[0 0 1]; X+[0 1 0]; X+[0 1 1]; X+[1 0 0]; X+[1 0 1]; X+[1 1 0]; X+[1 1 1]];
X = a*X;
Z = 4*ones(1, 64);

Ns=16*8; %# Number of states

setup(X, Ns, Z, S, R, true);
end