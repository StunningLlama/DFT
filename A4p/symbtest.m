syms X1 A x0 b c Q n
X0 = x0*x0+b*x0+c;
rs=(4*pi/3*n).^(-1/3); %# Added internal conversion to rs
x=sqrt(rs); X=x.*x+b*x+c;
dx=(0.5)./x; %# Chain rule needs dx/drho!

out0=dx.*( ...
2*X1./(rs.*x)+A*( ...
(2)./x-(2*x+b)./X-4*b./(Q*Q+(2*x+b).*(2*x+b)) ...
-(b*x0)/X0*((2)./(x-x0)-(2*x+b)./X-4*(2*x0+b)./ ...
(Q*Q+(2*x+b).*(2*x+b)) ) ...
) ...
);
out=(-rs./(3*n)).*out0;

out2=-X1./rs ...
+ A*( ...
+log(x.*x./X)+2*b/Q*atan(Q./(2*x+b)) ...
-(b*x0)/X0*( ...
log((x-x0).*(x-x0)./X)+2*(2*x0+b)/Q*atan(Q./(2*x+b)) ...
) ...
);

outd = diff(out,n)
out2d = diff(out2,n)

simplify(out2d-out, 'IgnoreAnalyticConstraints', true)