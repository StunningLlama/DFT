function hw3()
syms x r s a f(x,r) f2(x,r,s) f3(x,r,a);
f = r^2-x^3-r*atan(x);
bifurcationdiagram(f, [-0.3 0.3 -0.3 0.3], r, "4210_test", "")
f2 = s-r*x+x^2/(1+x^2);
f2tmp = subs(f2, r, 0.6);
f2tmp2 = subs(f2, r, 0.4);
bifurcationdiagram(f2tmp, [-0.5 1.5 -0.5 4], s, "4210_3_1", ", r=0.6")
bifurcationdiagram(f2tmp2, [-0.5 1.5 -0.5 4], s, "4210_3_2", ", r=0.4")
stabilitydiagram(f2, [-1 1 -1 1], r, s, [-1 3 1/20], "4210_3_3", "");

%f3 = r*x-a*x^2-x^3;
%stabilitydiagram(f3, [-3 3 -3 3], r, a, [-3 1 1/20]);

p4plots1();
p4plots2();
end

function saveimg(name)
saveas(gcf, "C:\Users\Brandon\Desktop\2218dia\" + name);
end

function bifurcationdiagram(func, limits, parameter, name, titlename)
syms x r;
rmin = limits(1);
rmax = limits(2);
xmin = limits(3);
xmax = limits(4);
f = matlabFunction(func, 'Vars', [x, parameter]);

findnumberofroots(@(x,r)(f(x,0.0948)), [0:0.1:2]);
g = matlabFunction(diff(func,x), 'Vars', [x, parameter]);
[R X] = meshgrid([0:0.002:1]*(rmax-rmin)+rmin,[0:0.002:1]*(xmax-xmin)+xmin);
[C h] = contourf(R, X, sign(g(X,R))+0.001*sin(X), '');
caxis([-2 2]);
set(h,'LineColor','none');
colormap(cool)
hold on;
fimplicit(@(x,y)(f(y,x)), [rmin rmax xmin xmax], 'white');
hold off;
xlabel(string(parameter));
ylabel("x");
title("Bifurcation diagram, Blue = stable, purple = unstable" + titlename);
saveimg(name + ".png");
end

function stabilitydiagram(func, limits, parameterx, parametery, searchparams, name, titlename)
syms x;
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
f = matlabFunction(func, 'Vars', [parameterx, parametery, x]);
[X Y] = meshgrid([0:0.002:1]*(xmax-xmin)+xmin,[0:0.002:1]*(ymax-ymin)+ymin);
ROOTS = X*0.0;

lowexponent = searchparams(1);
highexponent = searchparams(2);
division = searchparams(3);

for yi = 1:size(X, 1)
    for xi = 1:size(X,2)
        ROOTS(yi,xi) = findnumberofroots(@(x) f(X(yi,xi), Y(yi,xi),x), [-(10.^[highexponent:-division:lowexponent]) 0 10.^[lowexponent:division:highexponent]]);
    end
end

[C h] = contourf(X, Y, ROOTS, '');
caxis([0 3]);
set(h,'LineColor','none');
colormap(copper)
colorbar;
xlabel(string(parameterx));
ylabel(string(parametery));
title("Stability diagram colored by number of fixed points" + titlename);
saveimg(name + ".png");
end

function n = findnumberofroots(f, range)
vec = sign(f(range));
%vec(vec==0) = 1;
s = size(vec, 2);
n = sum(-vec(1:s-1) .* vec(2:s) + 1)*0.5;
end

function p4plots1
limits = [-5 5 -5 5]
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
[X Y] = meshgrid([0:0.05:1]*(xmax-xmin)+xmin,[0:0.05:1]*(ymax-ymin)+ymin);
U = (-9*X + 15*Y)/4;
V = (-15*X + 9*Y)/4;
t = [-5:0.1:5];
xp1 = t;
yp1 = t;
xp2 = t;
yp2 = -t;
quiver(X,Y,U,V, 2);
hold on;
plot(xp1, yp1, xp2, yp2);
hold off;
title("Phase portrait")
xlabel("x");
ylabel("y");
saveimg("4210_3_4.png");
end


function p4plots2
limits = [-5 5 -5 5]
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
[X Y] = meshgrid([0:0.05:1]*(xmax-xmin)+xmin,[0:0.05:1]*(ymax-ymin)+ymin);
U = (0*X + 1*Y)/4;
V = (-4*X - 2*Y)/4;
t = [-5:0.1:5];
t2 = [-3.3:0.1:3.3];
xp1 = t;
yp1 = 0*t;
xp2 = t2;
yp2 = -1.5*t2;
quiver(X,Y,U,V, 2);
hold on;
plot(xp1, yp1, xp2, yp2);
hold off;
title("Phase portrait")
xlabel("x");
ylabel("y");
saveimg("4210_3_5.png");
end