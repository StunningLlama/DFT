%# Code to solve Poisson's equation
%# Compute distances dr to center point in cell
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
fprintf('Ewald energy: %20.16f\n',Unum-Uself);