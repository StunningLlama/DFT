%# Code to solve Poisson's equation
%# Compute distances dr to center point in cell
function E = ewald()
global gbl_S; global gbl_R; global gbl_G2; global gbl_X; global gbl_G; global gbl_Sf; global gbl_X; global gbl_f; global gbl_r; global gbl_Z;
S = gbl_S;
R = gbl_R;
G2 = gbl_G2;
X = gbl_X;
G = gbl_G;
Sf = gbl_Sf;
X = gbl_X;
r = gbl_r;
Z = gbl_Z;
dr= sqrt(sum((r - ones(prod(S), 1)*sum(R,2)'/2).^2, 2)); %# <=== CODE INSERTION # 1
%# Compute two normalized Gaussians (widths 0.50 and 0.75)
sigma1=0.25;
g1=Z*exp(-dr.^2/(2*sigma1^2))/sqrt(2*pi*sigma1^2)^3;

%# Use structure factor to create all atoms
n=cI(cJ(g1).*Sf); n=real(n);

%# Check norms and integral (should be near 1 and 0, respectively)
fprintf('Normalization check on g1: %20.16f\n',sum(g1)*det(R)/prod(S));
fprintf('Total charge check: %20.16f\n',sum(n)*det(R)/prod(S));

phi=cI(Linv(-4*pi*O(cJ(n)))); %# <=== CODE INSERTION # 2
phi=real(phi);

Unum=0.5*real(cJ(phi)'*O(cJ(n)));
Uself=Z^2/(2*sqrt(pi))*(1/sigma1)*size(X,1);
E = Unum - Uself;
end