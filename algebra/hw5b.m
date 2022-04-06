function main
problem1();
end

function saveimg(name)
saveas(gcf, "C:\Users\Brandon\Desktop\2218dia\" + name);
end

function problem1
syms xdot ydot x y mu;
format compact;
xdot = -mu*y+x*y;
ydot = mu*x+(x^2-y^2)/2;
sol = solve(xdot==0, ydot==0, x, y);
j11 = diff(xdot,x);
j12 = diff(xdot,y);
j21 = diff(ydot,x);
j22 = diff(ydot,y);
det = j11*j22-j12*j21;
trace = j11+j22;
%sol.x
%sol.y
for i = [1:4]
sol.x(i)
sol.y(i)
delta = simplify(subs(det, [x,y], [sol.x(i), sol.y(i)]))
tau = simplify(subs(trace, [x,y], [sol.x(i), sol.y(i)]))
difference = simplify(expand(tau^2-4*delta), 'IgnoreAnalyticConstraints', true)
end
end