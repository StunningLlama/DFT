%function relax()

X=[8 8 8; 8 8 8; 8 8 8; 8 8 8; 8 8 8] + [0 0 0; 0 1.5 -0.3; 1 -0.5 -0.3; -1 -0.5 -0.3; 0 1.5 0];
setup(X, 5, [6 1 1 1 1]);
numatoms = 5;

%X=[8 8 8; 8+1.55 8 8];
%setup(X, 1, 1);
%numatoms = 2;

W = iterate();

dW = [];
dX = [];
setupPccgWavefunc(W);

for a = 1:(3*numatoms)
    dXi = zeros(3*numatoms, 1);
    dXi(a)=1;
    dXi = reshape(dXi, [numatoms, 3]);
    dX(:,:,a) = dXi;
    dWi = pccgWavefunc(W, dXi, 50, 1);
    dW(:,:,a) = dWi;
end

K = zeros(3*numatoms, 3*numatoms);

for a = 1:(3*numatoms)
    for b = a:(3*numatoms)
        kab = calcSpringConstant(W, X, dX(:,:,a), dX(:,:,b), dW(:,:,a), dW(:,:,b))-getdsqEwald(X, dX(:,:,a), dX(:,:,b))
        K(a,b)=kab;
        K(b,a)=kab;
    end
end

F = reshape(getForces(W),[],1);
%X = X+reshape(inv(K)*F, [numatoms, 3]);
%end