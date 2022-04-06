function main
problem2();
end

function saveimg(name)
saveas(gcf, "C:\Users\Brandon\Desktop\2218dia\" + name);
end

function problem1
syms xdot ydot x y A1 A2 alpha;
xdot = x*(1-x)-A1*x*y;
ydot = alpha*y*(1-y)-A2*x*y;
sol = solve(xdot==0, ydot==0, x, y);
j11 = diff(xdot,x);
j12 = diff(xdot,y);
j21 = diff(ydot,x);
j22 = diff(ydot,y);
det = j11*j22-j12*j21;
trace = j11+j22;
for i = [1:4]
sol.x(i)
sol.y(i)
delta = simplify(subs(det, [x,y], [sol.x(i), sol.y(i)]))
tau = simplify(subs(trace, [x,y], [sol.x(i), sol.y(i)]))
simplify(expand(tau^2-4*delta), 'IgnoreAnalyticConstraints', true)
end
end

function problem2
figure(1);
limits = [-0.1 1.1 -0.1 1.1]
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
[X Y] = meshgrid([0:0.002:1]*(xmax-xmin)+xmin,[0:0.002:1]*(ymax-ymin)+ymin);
A1 = 1.2;
A2 = 1.2;
alpha = 3;
U = X.*(1-X)-A1*X.*Y;
V = alpha*Y.*(1-Y)-A2*X.*Y;
streamslice(X,Y,U,V);
a1 = (-1+A2-alpha)/A1
a2 = -a1*(1-A2+A1*a1-a1*alpha)/(2-A2+alpha+3*A1*a1);
hold on;
fplot(@(x)((x-1)*a1+a2*(x-1).^2), [0 1], 'r');
hold off;
axis(limits);
xlabel('x');
ylabel('y');
title('Competing species model, unstable manifold');
saveimg("4210_4_5.png");

figure(2);
limits = [-0.1 1.1 -0.1 1.1]
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
[X Y] = meshgrid([0:0.002:1]*(xmax-xmin)+xmin,[0:0.002:1]*(ymax-ymin)+ymin);
A1 = 0.7;
A2 = 0.7;
alpha = 1;
U = X.*(1-X)-A1*X.*Y;
V = alpha*Y.*(1-Y)-A2*X.*Y;
streamslice(X,Y,U,V);
axis(limits);
xlabel('x');
ylabel('y');
title('Competing species model, case (i)');
saveimg("4210_4_6.png");

figure(3);
limits = [-0.1 1.1 -0.1 1.1]
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
[X Y] = meshgrid([0:0.002:1]*(xmax-xmin)+xmin,[0:0.002:1]*(ymax-ymin)+ymin);
A1 = 0.7;
A2 = 1.4;
alpha = 1;
U = X.*(1-X)-A1*X.*Y;
V = alpha*Y.*(1-Y)-A2*X.*Y;
streamslice(X,Y,U,V);
axis(limits);
xlabel('x');
ylabel('y');
title('Competing species model, case (i)');
saveimg("4210_4_7.png");

figure(4);
limits = [-0.1 1.1 -0.1 1.1]
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
[X Y] = meshgrid([0:0.002:1]*(xmax-xmin)+xmin,[0:0.002:1]*(ymax-ymin)+ymin);
A1 = 1.4;
A2 = 0.7;
alpha = 1;
U = X.*(1-X)-A1*X.*Y;
V = alpha*Y.*(1-Y)-A2*X.*Y;
streamslice(X,Y,U,V);
axis(limits);
xlabel('x');
ylabel('y');
title('Competing species model, case (i)');
saveimg("4210_4_8.png");

figure(5);
limits = [-0.1 1.1 -0.1 1.1]
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
[X Y] = meshgrid([0:0.002:1]*(xmax-xmin)+xmin,[0:0.002:1]*(ymax-ymin)+ymin);
A1 = 1.4;
A2 = 1.4;
alpha = 1;
U = X.*(1-X)-A1*X.*Y;
V = alpha*Y.*(1-Y)-A2*X.*Y;
streamslice(X,Y,U,V);
axis(limits);
xlabel('x');
ylabel('y');
title('Competing species model, case (i)');
saveimg("4210_4_9.png");
end