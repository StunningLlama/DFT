mu=2
a=1.1
[T, Y] = ode45(@(t,x)([mu*(x(2)+x(1)-x(1).^3/3); (a-x(1))/mu]), [0 10], [1.09; -0.84]);
[T2, Y2] = ode45(@(t,x)([mu*(x(2)+x(1)-x(1).^3/3); (a-x(1))/mu]), [0 10], [1.09; -1.03]);
plot(T, Y(:,1),T2, Y2(:,1))
xlabel("Time");
ylabel("x");
title("x vs time for small and large disturbances")