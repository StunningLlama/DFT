syms a1 a2 u alpha A1 A2;
eqn = (a1*u + a2*u^2 + 1000.0*(u^3))*(alpha*(1-a1*u - a2*u^2 + 1000.0*(u^3)) - A2*(u+1)) +(u+1)*(u + A1*(a1*u + a2*u^2 + 1000.0*(u^3)))*(a1 + 2*a2*u + 1000.0*(u^2));
eq1 = collect(simplify(expand(eqn)), u);
eq1 = subs(eq1, u^6, 0);
eq1 = subs(eq1, u^5, 0);
eq1 = subs(eq1, u^4, 0);
eq1 = subs(eq1, u^3, 0);