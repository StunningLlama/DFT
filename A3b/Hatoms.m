global gbl_Vdual;
%# Set the orbital occupancies
global gbl_f;
global gbl_active;
gbl_f=1; %# The usual case of a constant filling of two electrons per orbital

Vps=-4*pi*Z./G2; Vps(1)=0.;
Vdual=cJ(Vps.*Sf);
gbl_Vdual = Vdual;

%# Finite difference test
Ns=1; %# Number of states
randn('seed',0.2004);
W=(randn(length(gbl_active),Ns)+i*randn(length(gbl_active),Ns));
%more off; %# View output as it is computed
%fdtest(W);

format long
%# Converge
W=sd(W,20); %# 20 iterations of simple sd() to get nearer to the minimum
W=W*inv(sqrtm(W'*O(W))); %# Restart as orthonormal functions
[W, Elist]=pccg(W,50,1); %# 50 iterations of pclm from same W
Elist(50)

[Psi, epsilon]=getPsi(W);
for st=1:Ns
fprintf('=== State # %d, Energy = %f ===\n',st,epsilon(st));
dat=abs(cI(Psi(:,st))).^2;
tmpS = [9; 1; 1];
for k=1:3
sl=slice(dat,S,tmpS(k),k);
name=sprintf("psi%d_m%d.ppm",st,k);
%ppm(name,sl*0.3,sl,sl);% system(["display " name "&"]);
contourf(sl);
pause();
end
end