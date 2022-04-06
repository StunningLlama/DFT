syms x y f(x,y)

f(x,y) = x/(1+x^2+y^2)

xdot = diff(f,x);
ydot = diff(f,y);
sol = solve(xdot==0, ydot==0, x, y);
j11 = diff(xdot,x);
j12 = diff(xdot,y);
j21 = diff(ydot,x);
j22 = diff(ydot,y);
d2 = [j11 j12; j21 j22]
p1 = subs(d2, [x y], [1 0]);
p2 = subs(d2, [x y], [-1 0]);