h=0.00001;
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[T, Y] = ode45(@(t,x)([(-1+2*(x(1).^2+x(2).^2)-(x(1).^2+x(2).^2).^(3/2)).*x(1)-2*pi*x(2); ...
    (-1+2*(x(1).^2+x(2).^2)-(x(1).^2+x(2).^2).^(3/2)).*x(2)+2*pi*x(1)]), [0 1], [1+h, 0], options);
c1 = (Y(end,:)-[1 0])/h;

[T, Y2] = ode45(@(t,x)([(-1+2*(x(1).^2+x(2).^2)-(x(1).^2+x(2).^2).^(3/2)).*x(1)-2*pi*x(2); ...
    (-1+2*(x(1).^2+x(2).^2)-(x(1).^2+x(2).^2).^(3/2)).*x(2)+2*pi*x(1)]), [0 1], [1, h], options);

c2 = (Y2(end,:)-[1 0])/h;
mat = [c1' c2']
% 
% figure(1);
% limits = [-2 2 -2 2]
% xmin = limits(1);
% xmax = limits(2);
% ymin = limits(3);
% ymax = limits(4);
% [X Y] = meshgrid([0:0.002:1]*(xmax-xmin)+xmin,[0:0.002:1]*(ymax-ymin)+ymin);
% U = (-1+(X.^2+Y.^2)-(X.^2+Y.^2).^(3/2)).*X-2*pi*Y;
% V = (-1+(X.^2+Y.^2)-(X.^2+Y.^2).^(3/2)).*Y+2*pi*X;
% streamslice(X,Y,U,V);