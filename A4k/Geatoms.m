global gbl_Vdual;
%# Set the orbital occupancies
global gbl_f;
global gbl_active;
Ns=4; %# Number of states
gbl_f=[2;2/3;2/3;2/3]; %# The usual case of a constant filling of two electrons per orbital

%# Ge pseudopotential
Z=4;
lambda=18.5;
rc=1.052;
Gm=sqrt(G2);
Vps=-2*pi*exp(-pi*Gm/lambda).*cos(Gm*rc).*(Gm/lambda)./(1-exp(-2*pi*Gm/lambda));
for n=0:4
Vps=Vps+(-1)^n*exp(-lambda*rc*n)./(1+(n*lambda./Gm).^2);
end
Vps=Vps.*4*pi*Z./Gm.^2*(1+exp(-lambda*rc))-4*pi*Z./Gm.^2;
n=[1:4];
Vps(1)=4*pi*Z*(1+exp(-lambda*rc))*(rc^2/2+1/lambda^2* ...
(pi^2/6+sum((-1).^n.*exp(-lambda*rc*n)./n.^2)));


Vdual=cJ(Vps.*Sf);
gbl_Vdual = Vdual;

%# Finite difference test
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
%# Compute density for final viewing
Y=W*inv(sqrtm(W'*O(W))); %# Orthonormal wave functions
n=getn(cI(Y),gbl_f); %# Charge density
%# Basic slice from ‘‘100’’ edge of cell
sl=reshape(n(1:S(1)*S(2)),S(1),S(2));
%# Make and view image
contourf(sl);
%# 110 slice cutting bonds (assumes cube)
sl=reshape(n(find(M(:,2)==M(:,3))),S(1),S(2));
%# Expand by 2, drop data to restore (approximate) aspect ratio
li=find(rem([1:size(sl,1)],3)~=0); sl=sl(li,:);
%# Make and view image
figure(2);
contourf(sl);