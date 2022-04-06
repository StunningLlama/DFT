syms W Z r1 r2 r12 psi
psi = (1+W*r12)*exp(-Z*(r1+r2));
Lpsi = exp(-Z*(r1+r2))*((-2*W/(r12) + W*Z/(2*r12)*((r1^2-r2^2+r12^2)/r1+(r2^2-r1^2+r12^2)/r2) + (1+W*r12)*(Z/r1+Z/r2-Z^2)));
integrand = psi*(Lpsi-psi*(1/r1+1/r2-1/r12));
int1 = int(r12*integrand, r12, r2-r1, r2+r1)
int1 = simplify(int1, 'IgnoreAnalyticConstraints', true, 'Steps', 50)
int2 = int(r2*int1, r2, r1, inf)
int2 = simplify(int2, 'IgnoreAnalyticConstraints', true, 'Steps', 50)
int3 = int(r1*int2, r1, 0, inf)
int3 = simplify(int3, 'IgnoreAnalyticConstraints', true, 'Steps', 50)
int3