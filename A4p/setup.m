function setup(X, Ns, Z, S, R, kS, pseudopotential, nonlocal)
global gbl_S; global gbl_R; global gbl_G2; global gbl_Ns;
global gbl_G2c; global gbl_active; global gbl_X; global gbl_G; global gbl_Gc; global gbl_Sf; global gbl_M; global gbl_r; global gbl_Z;
global gbl_Vps;

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


global gbl_kvectors; global gbl_kpoints; global gbl_weights; global gbl_kS;

gbl_kS = kS;
kms=[0:prod(kS)-1]';
km1=rem(kms,kS(1));
km2=rem(floor(kms/kS(1)),kS(2));
km3=rem(floor(kms/(kS(1)*kS(2))),kS(3));
kM=[km1, km2, km3];
kR = kM*inv(diag(kS));

gbl_kpoints = prod(kS);
gbl_kvectors = kR*2*pi*inv(R);
gbl_weights = ones(1,gbl_kpoints)/gbl_kpoints;

%# Locate edges (assume S’s are even!) and determine max ’ok’ G2
if any(rem(S,2)~=0)
    fprintf("Odd dimension in S, cannot continue...\n");
    return;
end
eS=S/2+0.5;
edges=find(any(abs(M-ones(size(M,1),1)*eS')<1,2));

%# Compute active list and corresponding G2’s
G2mx=min(G2(edges));
gbl_Ecutoff = G2mx/4;

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

global gbl_gperm;
gbl_gperm = [1:prod(S)]';

%# Computation of structure factor
chargefactor = ones(size(G, 1), 1)*Z;
Sf=sum(chargefactor.*exp(-i*G*X'), 2);

gbl_S=S; gbl_R=R; gbl_G2=G2; gbl_G = G; gbl_Sf = Sf; gbl_r = r;


global gbl_Vdual;
%# Set the orbital occupancies
global gbl_f;
%gbl_f=1; %# The usual case of a constant filling of two electrons per orbital

gbl_f=2*ones(Ns,1);
gbl_Ns = Ns;

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
    
    gbl_Vps = vps;
    Vdual=cJdag(gbl_Vps.*Sf);
    gbl_Vdual = Vdual;
else
    gbl_Vps=-4*pi./G2; gbl_Vps(1)=0.;
    Vdual=cJdag(gbl_Vps.*Sf);
    gbl_Vdual = Vdual;
end

global gbl_K;
global gbl_Vnl;
global gbl_Mnl;

if (nonlocal)
    gbl_Vps = 0*gbl_Vps;
    gbl_Vdual = 0*gbl_Vdual;
    gbl_Mnl = diag([-1]);
    dr= sqrt(sum((r - ones(prod(S), 1)*sum(R,2)'/2).^2, 2)); %# <=== CODE INSERTION # 1
    %# Compute two normalized Gaussians (widths 0.50 and 0.75)
    sigma1=0.25;
    vtmp = exp(-dr.^2/(2*sigma1^2))/sqrt(2*pi*sigma1^2)^3;
    Sf3=sum(exp(-i*gbl_G*(-sum(R,2)/2)), 2);
    gbl_Vnl = real(cI(cJ(vtmp).*Sf3,1));
    %gbl_Vnl = vtmp;
    chargefactor = ones(size(gbl_Gc{1}, 1), 1)*Z;
    Sf2=sum(chargefactor.*exp(-i*gbl_Gc{1}*X'), 2);
    gbl_K = O(cJcomp(gbl_Vnl, 1).*Sf2);
end
end