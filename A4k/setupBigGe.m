function setupBigGe(resolution)

S=[resolution*2; resolution*2; resolution*2];
a=5.66/0.52917721; %# Lattice constant (converted from angstroms to bohrs)
R=2*a*diag(ones(3,1));

X=[0.00 0.00 0.00 %# diamond lattice in cubic cell
0.25 0.25 0.25
0.00 0.50 0.50
0.25 0.75 0.75
0.50 0.00 0.50
0.75 0.25 0.75
0.50 0.50 0.00
0.75 0.75 0.25];

X = [X; X+[0 0 1]; X+[0 1 0]; X+[0 1 1]; X+[1 0 0]; X+[1 0 1]; X+[1 1 0]; X+[1 1 1]];
X = a*X;
Z = 4*ones(1, 64);

global gbl_S; global gbl_R; global gbl_G2; global gbl_Ns;
global gbl_G2c; global gbl_active; global gbl_X; global gbl_G; global gbl_Gc; global gbl_Sf; global gbl_M; global gbl_r; global gbl_Z;
global gbl_Vps;
global gbl_Ecutoff;

gbl_X = X;
gbl_Z = Z;
ms=[0:prod(S)-1]';
m1=rem(ms,S(1));
m2=rem(floor(ms/S(1)),S(2));
m3=rem(floor(ms/(S(1)*S(2))),S(3));
M=[m1, m2, m3];
gbl_M = M;

n1=m1-(m1>S(1)/2)*S(1);
n2=m2-(m2>S(2)/2)*S(2);
n3=m3-(m3>S(3)/2)*S(3);
N=[n1,n2,n3];

r = M*inv(diag(S))*(R');
G = 2*pi*N*inv(R);
G2 = sum(G.^2, 2);

global gbl_kvectors; global gbl_kpoints; global gbl_weights;
gbl_kpoints = 1;
gbl_kvectors = [0 0 0]*2*pi*inv(R);
gbl_weights = [1];


gbl_G2c = {};
gbl_active = {};
gbl_Gc = {};
for k = [1:gbl_kpoints]
    kvec = gbl_kvectors(k,:);
    karray = ones(size(G,1),1)*kvec;
    active=find(sum((G+karray).^2, 2)<gbl_Ecutoff); %# Sphere is 1/2 size (but looking at G^2!)
    G2c=G2(active);
    fprintf("Compression: %f (theoretical: %f)\n", ...
        length(G2)/length(G2c), 1/(4*pi*(1/4)^3/3));
    gbl_G2c{k} = G2c;
    gbl_active{k} = active;
    gbl_Gc{k} = G(active,:);
end

%# Computation of structure factor
chargefactor = ones(size(G, 1), 1)*Z;
Sf=sum(chargefactor.*exp(-i*G*X'), 2);

gbl_S=S; gbl_R=R; gbl_G2=G2; gbl_G = G; gbl_Sf = Sf; gbl_r = r;


global gbl_Vdual;
%# Set the orbital occupancies
global gbl_f;

% gbl_kpoints = 27;
% gbl_kvectors = K*2*pi*inv(R);
% gbl_weights = ones(1, 27)/27;

Ns=16*8; %# Number of states
gbl_Ns = Ns;
gbl_f=2*ones(16*8,1); %# The usual case of a constant filling of two electrons per orbital

lambda=18.5;
rc=1.052;
Gm=sqrt(G2);
Vps=-2*pi*exp(-pi*Gm/lambda).*cos(Gm*rc).*(Gm/lambda)./(1-exp(-2*pi*Gm/lambda));
for n=0:4
Vps=Vps+(-1)^n*exp(-lambda*rc*n)./(1+(n*lambda./Gm).^2);
end
Vps=Vps.*4*pi./Gm.^2*(1+exp(-lambda*rc))-4*pi./Gm.^2;
n=[1:4];
Vps(1)=4*pi*(1+exp(-lambda*rc))*(rc^2/2+1/lambda^2* ...
(pi^2/6+sum((-1).^n.*exp(-lambda*rc*n)./n.^2)));

Vdual=cJ(Vps.*Sf);
gbl_Vdual = Vdual;

end