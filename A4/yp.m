
function main
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[tout, yout] = ode45(@f, [0, 10], [0; 0], opts);
[tout2, yout2] = ode45(@f2, [0, 10], [0], opts);
%plot(tout, yout.^3.*(ones(size(tout))*[2 1]), tout2, yout2.^3, tout, inp(tout)>-69);
plot(tout, yout(:,1).^3*4, tout2, yout2.^4, tout, inp(tout)>-69);
axis([-inf, inf, 0, 1]);
end

function val = inp(t)
dt1 = 3;
dt2 = 6;
current = ((t < 1+dt1) & (t > 1) | (t < 1+dt1+dt2) & (t > 1+dt2));
val = -70*(1-current) + 0*current;
end

function y=sigmoid(x)
y = 1/(1+exp(-4*x));
end

function yp = f(t, y)
openrate = sigmoid((inp(t)+15)/80);
blockrate = sigmoid((inp(t)+50)/20);

rtc = 0.3;
iatc = 1.2;
diatc = 2.5;
btc = 0.6;

ra = openrate/rtc; %rising time
rb = (1-openrate)/btc;
rc = blockrate/iatc; %inactivation time
rd = (1-blockrate)/diatc; %deinactivation time
yp = [ra*(1-y(1)-y(2))-rb*y(1) + rd*y(2) -  rc*y(1); rc*y(1) - rd*y(2)];
end

function yp = f2(t, y)
openrate = sigmoid((inp(t)+15)/80);

atc = 0.8;
btc = 1.5;

ra = openrate/atc;
rb = (1-openrate)/btc;
yp = [ra*(1-y(1)) - rb*(y(1))];
end