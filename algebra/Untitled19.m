a=0.2;
b=0.2;
c=5.7;
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

ic = 5*[randn() randn() randn()]+[0 0 5];
[T, Y] = ode45(@(t,x)([-x(2)-x(3); x(1)+a*x(2); b+x(3).*(x(1)-c)]), [0 10000], ic, options);
plot3(Y(:,1), Y(:,2), Y(:,3));
pause();

for n = [1:3]
ic = 5*[randn() randn() randn()]+[0 0 5];
[T, Y] = ode45(@(t,x)([-x(2)-x(3); x(1)+a*x(2); b+x(3).*(x(1)-c)]), [0 100], ic, options);
plot3(Y(:,1), Y(:,2), Y(:,3));
hold on;
end
hold off;
pause();

ic = [-4.20 0.69 0.1];
icx = ic + [1e-4 0 0];
icy = ic + [0 1e-4 0];
icz = ic + [0 0 1e-4];
[T, Y0] = ode45(@(t,x)([-x(2)-x(3); x(1)+a*x(2); b+x(3).*(x(1)-c);...
    -x(5)-x(6); x(4)+a*x(5); b+x(6).*(x(4)-c);...
    -x(8)-x(9); x(7)+a*x(8); b+x(9).*(x(7)-c);...
    -x(11)-x(12); x(10)+a*x(11); b+x(12).*(x(10)-c)]...
), [0 500], [ic icx icy icz], options);
Y = Y0(:,1:3);
Yx = Y0(:,4:6);
Yy = Y0(:,7:9);
Yz = Y0(:,10:12);
dx = sqrt(sum((Yx-Y).^2, 2));
dy = sqrt(sum((Yy-Y).^2, 2));
dz = sqrt(sum((Yz-Y).^2, 2));

plot(T, log(dx));
hold on;
fplot(@(x)(0.075*x-9.5));
hold off;
title('Log distance vs time for x-displaced trajectory');
xlabel('t');
ylabel('Ln(dist)');
legend('', 'y=0.075x-9.5');

pause();
plot(T, log(dy));
hold on;
fplot(@(x)(0.075*x-12));
hold off;
title('Log distance vs time for y-displaced trajectory');
xlabel('t');
ylabel('Ln(dist)');
legend('', 'y=0.075x-12');

pause();
plot(T, log(dz));
hold on;
fplot(@(x)(0.075*x-12));
hold off;
title('Log distance vs time for z-displaced trajectory');
xlabel('t');
ylabel('Ln(dist)');
legend('', 'y=0.075x-12');


Y_x = Y0(:,1);
Yx_x = Y0(:,4);
Yy_x = Y0(:,7);
Yz_x = Y0(:,10);

tiledlayout(3,1)
nexttile
plot(T, Y_x);
hold on;
plot(T, Yx_x);
hold off;
title('Time series plot');
xlabel('t');
ylabel('x');
legend('reference', 'perturbed x');

nexttile
plot(T, Y_x);
hold on;
plot(T, Yy_x);
hold off;
title('Time series plot');
xlabel('t');
ylabel('x');
legend('reference', 'perturbed y');

nexttile
plot(T, Y_x);
hold on;
plot(T, Yz_x);
hold off;
title('Time series plot');
xlabel('t');
ylabel('x');
legend('reference', 'perturbed z');