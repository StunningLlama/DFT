%warning("off", "Octave:divide-by-zero");

%function relax()

%  bondlen = 108.7*0.0188972599;
%  angle0 = deg2rad(109.5);
%  angle1 = deg2rad(0);
% angle2 = deg2rad(120);
% angle3 = deg2rad(240);
% X=[8 8 8; 8 8 8; 8 8 8; 8 8 8; 8 8 8] + bondlen*[0 0 0; cos(angle1) sin(angle1) cos(angle0); cos(angle2) sin(angle2) cos(angle0); cos(angle3) sin(angle3) cos(angle0); 0 0 1];
% Z = [6 1 1 1 1];
% nstates = 5;
% setup(X, nstates, Z, true);
% numatoms = 5;

X=[8 8 8; 8+1.55 8 8];
nstates = 1;
Z = [1];
setup(X, nstates, Z, true);
numatoms = 2;

W = iterate(20);
visualize(W, X);
E0 = getE(W) + ewald();

tic
%Calculate spring const analytically
dW = [];
dX = [];
setupPccgWavefunc(W);

for a = 1:(3*numatoms)
    dXi = zeros(3*numatoms, 1);
    dXi(a)=1;
    dXi = reshape(dXi, [numatoms, 3]);
    dX(:,:,a) = dXi;
    dWi = pccgWavefunc(W, dXi, 50, 1);
    dW(:,:,:,a) = dWi;
end

Kanal = zeros(3*numatoms, 3*numatoms);

for a = 1:(3*numatoms)
    for b = a:(3*numatoms)
        kab = calcSpringConstant(W, X, dX(:,:,a), dX(:,:,b), dW(:,:,:,a), dW(:,:,:,b))-getdsqEwald(X, dX(:,:,a), dX(:,:,b))
        Kanal(a,b)=kab;
        Kanal(b,a)=kab;
    end
end
t1 = toc

tic
Wp = [];
E = [];
h=0.0001;
Knum = zeros(3*numatoms, 3*numatoms);

for a = 1:(3*numatoms)
    dXi = zeros(3*numatoms, 1);
    dXi(a)=1;
    dXi = reshape(dXi, [numatoms, 3]);
    setup(X+h*dXi, nstates, Z, true);
    Wp(:,:,:,a) = pccg(W,50,1);
    E(a) = getE(Wp(:,:,:,a)) + ewald();
end

for a = 1:(3*numatoms)
    for b = a:(3*numatoms)
        dXi = zeros(3*numatoms, 1);
        if (a == b)
            dXi(a) = 2;
        else
            dXi(a)=1;
            dXi(b)=1;
        end
        dXi = reshape(dXi, [numatoms, 3]);
        setup(X+h*dXi, nstates, Z, true);
        Wab = pccg(W,50,1);
        Eab = getE(Wab) + ewald();
        kab = -(Eab-E(a)-E(b)+E0)/h^2;
        Knum(a,b)=kab;
        Knum(b,a)=kab;
    end
end
t2 = toc

Kanal
Knum

F = reshape(getForces(W),[],1);
%X = X+reshape(inv(K)*F, [numatoms, 3]);
%end