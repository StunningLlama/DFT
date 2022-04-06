syms r theta psi1 psi2 E

psi1 = (2-r)*exp(-r/2)/((sym(32)*pi)^{1/2});
psi2 = r*exp(-r/2)*cos(theta)/((sym(32)int*pi)^{1/2});
psi1Hpsi2 = psi1*E*r*cos(theta)*psi2;
int1 = int(2*pi*r^2*sin(theta)*psi1Hpsi2, theta, 0, pi);
int2 = int(int1, r, 0, inf);
