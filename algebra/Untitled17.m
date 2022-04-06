%for rho = [15:0.5:50]
%for sigma = [5:0.5:15]
for j = [1]
    rho = 25;
    sigma = 10;
    beta = 8/3;
    ic = 10*[randn() randn() randn()];
    [T, Y] = ode45(@(t,x)([sigma*(x(2)-x(1)); x(1)*(rho-x(3))-x(2); x(1)*x(2)-beta*x(3)]), [0 10000], ic);
    plot3(Y(:,1), Y(:,2), Y(:,3),'Color',[1, 0, 0, 0.2]);
    axis([-20.0000   25.0000  -25.4223   29.0956         0   60.0000]);
    pause(0.2);
end
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