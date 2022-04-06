global gbl_S; global gbl_R; global gbl_G2;
global gbl_G2c; global gbl_active;

S = [32; 32; 32];

R=diag([8; 8; 8]);
%# Define atomic locations and nuclear charge
X=[0 0 0]; Z=1;

ms=[0:prod(S)-1]';
m1=rem(ms,S(1));
m2=rem(floor(ms/S(1)),S(2));
m3=rem(floor(ms/(S(1)*S(2))),S(3));
M=[m1, m2, m3];

n1=m1-(m1>S(1)/2)*S(1);
n2=m2-(m2>S(2)/2)*S(2);
n3=m3-(m3>S(3)/2)*S(3);
N=[n1,n2,n3];

r = M*inv(diag(S))*(R');
G = 2*pi*N*inv(R);
G2 = sum(G.^2, 2);

%# Locate edges (assume S’s are even!) and determine max ’ok’ G2
if any(rem(S,2)~=0)
fprintf("Odd dimension in S, cannot continue...\n");
return;
end
eS=S/2+0.5;
edges=find(any(abs(M-ones(size(M,1),1)*eS')<1,2));

%# Compute active list and corresponding G2’s
G2mx=min(G2(edges));
active=find(G2<G2mx/4); %# Sphere is 1/2 size (but looking at G^2!)
G2c=G2(active);
fprintf("Compression: %f (theoretical: %f)\n", ...
length(G2)/length(G2c), 1/(4*pi*(1/4)^3/3));
gbl_G2c = G2c;
gbl_active = active;

%# Computation of structure factor
Sf=sum( exp(-i*G*X'), 2);

gbl_S=S; gbl_R=R; gbl_G2=G2;