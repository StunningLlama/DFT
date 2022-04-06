function testgreens()
Q = 3;
omega=sqrt(1-1/(4*Q^2));
tau=2*Q;
xvec = [-20:0.02:20];
yvec = [];
for x=xvec
    yvec = [yvec u(x,omega,tau)];
end
h = 0.02;
yvecd = 0.5*([diff(yvec,1)/h 0]+[0 diff(yvec,1)/h]);
yvecdd = [0 diff(yvec,2)/(h^2) 0];

plot(xvec, yvec, xvec,yvecdd+yvecd/Q+yvec, xvec, f(xvec));
end

function g = G(t, omega, tau)
    g = (t>0).*(omega+1/(omega*tau^2)).*exp(-t/tau).*sin(omega*t);
end

function y = u(t, omega, tau)
    y = integral(@(s)(G(t-s, omega, tau).*f(s)),-10,10,'RelTol',1e-10,'AbsTol',1e-10);
end


function g = G2(t, omega2, tau)
    g = (t>0).*(omega+1/(omega*tau^2)).*exp(-t/tau).*sin(omega*t);
end

function y = u(t, omega, tau)
    y = integral(@(s)(G(t-s, omega, tau).*f(s)),-10,10,'RelTol',1e-10,'AbsTol',1e-10);
end


function y = f(t)
%y=exp(-t.^4/2);
y=sin(t+5*sin(t));
end
