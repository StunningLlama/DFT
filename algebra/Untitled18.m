syms phidot psidot phi psi a
phidot = a*sin(psi)-2*a*sin(phi) - a*sin(psi+phi)
psidot = a*sin(phi)-2*a*sin(psi) - a*sin(psi+phi)
sol = solve(phidot == 0, psidot == 0)
sol
df = [diff(phidot, phi, 1), diff(phidot, psi, 1); diff(psidot, phi, 1), diff(psidot, psi, 1)]
detf = det(df)
trf = trace(df)
for i = [1:size(sol.phi)]
    [sol.phi(i) sol.psi(i) subs(detf,[phi psi],[sol.phi(i) sol.psi(i)]) subs(trf,[phi psi],[sol.phi(i) sol.psi(i)])]
end