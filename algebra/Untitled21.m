syms x t A B L delta
L = (-A*sin(t)+B*cos(t))^2-(A*cos(t)+B*sin(t))^2;
sol = int(L, t, 0, delta);
sol