global gbl_Vdual;
%# Set the orbital occupancies
global gbl_f;
gbl_f=2; %# The usual case of a constant filling of two electrons per orbital

omega = 2;
V = 0.5*omega^2*sum((r - ones(prod(S), 1)*sum(R,2)'/2).^2, 2);
gbl_Vdual = cJdag(O(cJ(V)));

%# Finite difference test
Ns=4; %# Number of states
randn('seed',0.2004);
W=(randn(prod(S),Ns)+i*randn(prod(S),Ns));
more off; %# View output as it is computed
fdtest(W);
%pause();

W = W*sqrtm(inv(W'*O(W)));
%W = sd(W, 250);

%# Allow for more digits in printouts
format long
%# Converge using steepest descents
W=sd(W,400);
%# Extract and display final results
[Psi, epsilon]=getPsi(W);
for st=1:Ns
fprintf('=== State # %d, Energy = %f ===\n',st,epsilon(st));
dat=abs(cI(Psi(:,st))).^2;
for k=1:3
sl=slice(dat,S,S(k)/2,k);
name=sprintf("psi%d_m%d.ppm",st,k);
%ppm(name,sl*0.3,sl,sl);% system(["display " name "&"]);
contourf(sl);
pause();
end
end