n = 3.34343;

X1 = 0.75*(3.0/(2.0*pi))^(2.0/3.0);
A = 0.0310907;
x0 = -0.10498;
b = 3.72744;
c = 12.9352;
Q = sqrt(4*c-b*b);
X0 = x0*x0+b*x0+c;
rs=(4*pi/3*n).^(-1/3); %# Added internal conversion to rs
x=sqrt(rs); X=x.*x+b*x+c;
excpp = (2.*X1.*((4.*pi.*n)./3).^(1./2) - A.*((4.*b)./(Q.^2 + (b + 2./((4.*n.*pi)./3).^(1./6)).^2) - 2.*((4.*pi.*n)./3).^(1./6) + (b + 2./((4.*n.*pi)./3).^(1./6))./(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3)) - (b.*x0.*(2./(x0 - 1./((4.*n.*pi)./3).^(1./6)) + (4.*b + 8.*x0)./(Q.^2 + (b + 2./((4.*n.*pi)./3).^(1./6)).^2) + (b + 2./((4.*n.*pi)./3).^(1./6))./(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3))))./(x0.^2 + b.*x0 + c)))./(6.*n.^2.*((4.*n.*pi)./3).^(1./6)) + (A.*((((4.*pi)./(9.*((4.*n.*pi)./3).^(4./3)) + (2.*b.*pi)./(9.*((4.*n.*pi)./3).^(7./6))).*(b + 2./((4.*n.*pi)./3).^(1./6)))./(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*n.*pi)./3).^(1./3)).^2 - (4.*pi)./(9.*((4.*n.*pi)./3).^(7./6).*(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3))) - (4.*pi)./(9.*((4.*n.*pi)./3).^(5./6)) + (b.*x0.*((4.*pi)./(9.*((4.*n.*pi)./3).^(7./6).*(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3))) + (4.*pi)./(9.*(x0 - 1./((4.*n.*pi)./3).^(1./6)).^2.*((4.*n.*pi)./3).^(7./6)) - (((4.*pi)./(9.*((4.*n.*pi)./3).^(4./3)) + (2.*b.*pi)./(9.*((4.*n.*pi)./3).^(7./6))).*(b + 2./((4.*n.*pi)./3).^(1./6)))./(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*n.*pi)./3).^(1./3)).^2 - (8.*pi.*(b + 2./((4.*n.*pi)./3).^(1./6)).*(4.*b + 8.*x0))./(9.*(Q.^2 + (b + 2./((4.*n.*pi)./3).^(1./6)).^2).^2.*((4.*n.*pi)./3).^(7./6))))./(x0.^2 + b.*x0 + c) + (32.*b.*pi.*(b + 2./((4.*n.*pi)./3).^(1./6)))./(9.*(Q.^2 + (b + 2./((4.*n.*pi)./3).^(1./6)).^2).^2.*((4.*n.*pi)./3).^(7./6))) - (4.*X1.*pi)./(3.*((4.*n.*pi)./3).^(1./2)))./(6.*n.*((4.*n.*pi)./3).^(1./6)) + (pi.*(2.*X1.*((4.*pi.*n)./3).^(1./2) - A.*((4.*b)./(Q.^2 + (b + 2./((4.*n.*pi)./3).^(1./6)).^2) - 2.*((4.*pi.*n)./3).^(1./6) + (b + 2./((4.*n.*pi)./3).^(1./6))./(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3)) - (b.*x0.*(2./(x0 - 1./((4.*n.*pi)./3).^(1./6)) + (4.*b + 8.*x0)./(Q.^2 + (b + 2./((4.*n.*pi)./3).^(1./6)).^2) + (b + 2./((4.*n.*pi)./3).^(1./6))./(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3))))./(x0.^2 + b.*x0 + c))))./(27.*n.*((4.*n.*pi)./3).^(7./6))

excp = - A.*(((4.*pi)./(9.*((4.*n.*pi)./3).^(4./3).*(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3))) - ((4.*pi)./(9.*((4.*n.*pi)./3).^(4./3)) + (2.*b.*pi)./(9.*((4.*n.*pi)./3).^(7./6)))./(((4.*n.*pi)./3).^(1./3).*(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*n.*pi)./3).^(1./3)).^2)).*((4.*pi.*n)./3).^(1./3).*(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3)) + (b.*x0.*((((((4.*pi)./(9.*((4.*n.*pi)./3).^(4./3)) + (2.*b.*pi)./(9.*((4.*n.*pi)./3).^(7./6))).*(x0 - 1./((4.*n.*pi)./3).^(1./6)).^2)./(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*n.*pi)./3).^(1./3)).^2 + (4.*pi.*(x0 - 1./((4.*n.*pi)./3).^(1./6)))./(9.*((4.*n.*pi)./3).^(7./6).*(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3)))).*(c + b./((4.*n.*pi)./3).^(1./6) + 1./((4.*pi.*n)./3).^(1./3)))./(x0 - 1./((4.*n.*pi)./3).^(1./6)).^2 + (4.*pi.*(2.*b + 4.*x0))./(9.*(Q.^2./(b + 2./((4.*n.*pi)./3).^(1./6)).^2 + 1).*(b + 2./((4.*n.*pi)./3).^(1./6)).^2.*((4.*n.*pi)./3).^(7./6))))./(x0.^2 + b.*x0 + c) - (8.*b.*pi)./(9.*(Q.^2./(b + 2./((4.*n.*pi)./3).^(1./6)).^2 + 1).*(b + 2./((4.*n.*pi)./3).^(1./6)).^2.*((4.*n.*pi)./3).^(7./6))) - (4.*X1.*pi)./(9.*((4.*n.*pi)./3).^(2./3))

excpp
excp

excpVWN(n)
excVWN(n)