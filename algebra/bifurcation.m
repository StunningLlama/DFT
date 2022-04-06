function main()
problem4a();
end

function saveimg(name)
saveas(gcf, "C:\Users\Brandon\Desktop\" + name);
end

function problem1()
rmin = -3;
rmax = 3;
xmin = -3;
xmax = 3;
[R X] = meshgrid([rmin:0.02:rmax],[xmin:0.02:xmax]);
[C h] = contourf(R, X, sign(g(X,R)), '');
caxis([-2 2]);
set(h,'LineColor','none');
colormap(cool)
hold on;
fimplicit(@(x,y)(f(y,x)), [rmin rmax xmin xmax], 'white');
hold off;
xlabel("r");
ylabel("x");
title("Bifurcation diagram, Blue = stable, purple = unstable");
saveimg("4210_2_3.png");
end


function problem1a()
syms x r f(x,r)
f = -x^3+2*x^2+r*x;
g = diff(f,x);
sol = solve(f==0, g==0);
[sol.r sol.x]
end

function xdot = f(x,r)
xdot = -x.^3+2*x.^2-r.*x;
end

function g = g(x,r)
g = -3*x.^2 + 4.*x - r;
end


function f2 = f2(x,r)
f2 = r.*x-sin(x);
end

function g2 = g2(x,r)
g2 = r-cos(x);
end

function problem2()
rmin = -3;
rmax = 3;
xmin = -10;
xmax = 10;
[R X] = meshgrid([rmin:0.02:rmax],[xmin:0.02:xmax]);
[C h] = contourf(R, X, sign(g2(X,R)), '');
caxis([-2 2]);
set(h,'LineColor','none');
colormap(cool)
hold on;
fimplicit(@(r,x)(f2(x,r)), [rmin rmax xmin xmax], 'white');
hold off;
xlabel("r");
ylabel("x");
title("Bifurcation diagram, Blue = stable, purple = unstable");
saveimg("4210_2_4.png");
end


function xdot = h(x, a, r)
    xdot = r.*x-a*x.^2-x.^3;
end

function xdot = g3(x, a, r)
    xdot = r-2*a*x-3*x.^2;
end

function problem3a()
rmin = -3;
rmax = 3;
xmin = -3;
xmax = 3;
a = -1;
[R X] = meshgrid([rmin:0.02:rmax],[xmin:0.02:xmax]);
[C h2] = contourf(R, X, sign(g3(X,a, R)), '');
caxis([-2 2]);
set(h2,'LineColor','none');
colormap(cool)
hold on;
fimplicit(@(r,x)(h(x, a, r)), [rmin rmax xmin xmax], 'white');
hold off;
xlabel("r");
ylabel("x");
title("Bifurcation diagram, Blue = stable, purple = unstable, a = " + num2str(a));
saveimg("4210_2_7.png");
end

function problem3()
    fimplicit3(@(x,y,z)h(z,x,y));%,'EdgeColor','none','FaceAlpha',.2
end


function problem4a()
syms x r f(x,r)
f = r*x-sin(x);
g = diff(f,x);
sol = solve(f==0, g==0);
[sol.r sol.x]
end
