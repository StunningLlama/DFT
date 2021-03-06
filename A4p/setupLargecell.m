function setupLargecell(X, Ns_small, Z, S_small, R_small, kS, pseudopotential)
global gbl_S; global gbl_R; global gbl_G2; global gbl_Ns;
global gbl_G2c; global gbl_active; global gbl_X; global gbl_G; global gbl_Gc; global gbl_Sf; global gbl_M; global gbl_r; global gbl_Z;
global gbl_Vps;


global gbl_kvectors; global gbl_kpoints; global gbl_weights; global gbl_kS;
global gbl_Ns_small;

gbl_Ns_small = Ns_small;

Ns_BIG = Ns_small*prod(kS);
S_BIG = S_small.*kS;
R_BIG = R_small*diag(kS);

gbl_kS = kS;
kms=[0:prod(kS)-1]';
km1=rem(kms,kS(1));
km2=rem(floor(kms/kS(1)),kS(2));
km3=rem(floor(kms/(kS(1)*kS(2))),kS(3));
kM=[km1, km2, km3];
kR = kM*inv(diag(kS));

gbl_kpoints = prod(kS);
gbl_kvectors = kR*2*pi*inv(R_small);
gbl_weights = ones(1,gbl_kpoints);%/gbl_kpoints;

ms_small=[0:prod(S_small)-1]';
m1_small=rem(ms_small,S_small(1));
m2_small=rem(floor(ms_small/S_small(1)),S_small(2));
m3_small=rem(floor(ms_small/(S_small(1)*S_small(2))),S_small(3));
M_small=[m1_small, m2_small, m3_small];
n1_small=m1_small-(m1_small>S_small(1)/2)*S_small(1);
n2_small=m2_small-(m2_small>S_small(2)/2)*S_small(2);
n3_small=m3_small-(m3_small>S_small(3)/2)*S_small(3);
N_small=[n1_small,n2_small,n3_small];

ms_BIG=[0:prod(S_BIG)-1]';
m1_BIG=rem(ms_BIG,S_BIG(1));
m2_BIG=rem(floor(ms_BIG/S_BIG(1)),S_BIG(2));
m3_BIG=rem(floor(ms_BIG/(S_BIG(1)*S_BIG(2))),S_BIG(3));
M_BIG=[m1_BIG, m2_BIG, m3_BIG];
n1_BIG=m1_BIG-(m1_BIG>S_BIG(1)/2)*S_BIG(1);
n2_BIG=m2_BIG-(m2_BIG>S_BIG(2)/2)*S_BIG(2);
n3_BIG=m3_BIG-(m3_BIG>S_BIG(3)/2)*S_BIG(3);
N_BIG=[n1_BIG,n2_BIG,n3_BIG];
gbl_M = M_BIG;

%r = M_small*inv(diag(S_small))*(R_small');
G0 = 2*pi*N_small*inv(R_small);
G20 = sum(G0.^2, 2);

%# Locate edges (assume S’s are even!) and determine max ’ok’ G2
if any(rem(S_small,2)~=0)
    fprintf("Odd dimension in S, cannot continue...\n");
    return;
end
eS=S_small/2+0.5;
edges=find(any(abs(M_small-ones(size(M_small,1),1)*eS')<1,2));

%# Compute active list and corresponding G2’s
G2mx=min(G20(edges));
gbl_Ecutoff = G2mx/4;

G2csmall = {};
activesmall = {};
Gcsmall = {};
Gfull = [];
G2full = [];
Coordfull = [];

totallength = zeros(1,gbl_kpoints);
activelength = zeros(1,gbl_kpoints);
totalindex = zeros(1,gbl_kpoints);
activeindex = zeros(1,gbl_kpoints);

global gbl_gperm; global gbl_activelength; global gbl_activeindex;
gbl_gperm = [];

for k = [1:gbl_kpoints]
    kvec = gbl_kvectors(k,:);
    karray = ones(size(G0,1),1)*kvec;
    G = G0+karray;
    G2 = sum(G.^2, 2);
    Coord = N_small*diag(kS) + kM(k,:);
    Gfull = [Gfull; G];
    G2full = [G2full; G2];
    Transformedcoord = Coord + (Coord<0)*diag(S_BIG);
    coordindices = coordtoindex(Transformedcoord, S_BIG)+1;
    gbl_gperm = [gbl_gperm; coordindices];
    
    active=find(G2<gbl_Ecutoff); %# Sphere is 1/2 size (but looking at G^2!)
    G2c=G2(active);
    fprintf("Compression: %f (theoretical: %f)\n", ...
        length(G2)/length(G2c), 1/(4*pi*(1/4)^3/3));
    totallength(k) = length(G2);
    activelength(k) = length(G2c);
    totalindex(k) = sum(totallength(1:k)) - totallength(k) + 1;
    activeindex(k) = sum(activelength(1:k)) - activelength(k) + 1;
    
    G2csmall{k} = G2c;
    activesmall{k} = coordindices(active);
    Gcsmall{k} = G(active,:);
end

G = Gfull;
G2 = G2full;
gbl_activelength = activelength;
gbl_activeindex = activeindex;

gbl_G2c = {};
gbl_active = {};
gbl_Gc = {};
for k = [1:gbl_kpoints]
    gbl_G2c{1}([activeindex(k): activeindex(k)+activelength(k)-1], 1) = G2csmall{k};
    %gbl_active{1}(activeindex(k) + [1:activelength(k)] - 1, 1) = activesmall{k}+totalindex(k)-1;
    gbl_active{1}([activeindex(k): activeindex(k)+activelength(k)-1], 1) = activesmall{k}
    gbl_Gc{1}([activeindex(k): activeindex(k)+activelength(k)-1], 1:3) = Gcsmall{k};
end

G_BIG = 2*pi*N_BIG*inv(R_BIG);
G2_BIG = sum(G_BIG.^2, 2);
AA = G_BIG(gbl_active{1},:);
BB = gbl_Gc{1};

gbl_r = M_BIG*inv(diag(S_BIG))*(R_BIG');

gbl_kS = [1;1;1];
kms=[0:prod(kS)-1]';
km1=rem(kms,kS(1));
km2=rem(floor(kms/kS(1)),kS(2));
km3=rem(floor(kms/(kS(1)*kS(2))),kS(3));
kM=[km1, km2, km3];
kR = kM*inv(diag(kS));

Xsmall = X;
Zsmall = Z;
X = [];
Z = [];
for k = [1:prod(kS)]
    offsets = kM*(R_small');
    X = [X; Xsmall + offsets(k,:)];
    Z = [Z Zsmall];
end
%TODO

gbl_kpoints = 1;
gbl_kvectors = [0 0 0];
gbl_weights = [1];

gbl_X = X;
gbl_Z = Z;

%TEST
% G = G_BIG;
% G2 = G2_BIG;
% gbl_G = G_BIG;
% gbl_G2 = G2_BIG;
% gbl_gperm = [1:prod(S_BIG)]';


%# Computation of structure factor
chargefactor = ones(size(G, 1), 1)*Z;
Sf=sum(chargefactor.*exp(-i*G*X'), 2);

gbl_S=S_BIG; gbl_R=R_BIG; gbl_G2=G2; gbl_G = G; gbl_Sf = Sf;



global gbl_Vdual;
%# Set the orbital occupancies
global gbl_f;
%gbl_f=1; %# The usual case of a constant filling of two electrons per orbital

gbl_f=2*ones(Ns_BIG,1);
gbl_Ns = Ns_BIG;

if (pseudopotential)
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
    
    Vdual=cJdag(Vps.*Sf);
    gbl_Vdual = Vdual;
else
    gbl_Vps=-4*pi./G2; gbl_Vps(1)=0.;
    Vdual=cJdag(gbl_Vps.*Sf);
    gbl_Vdual = Vdual;
end

end