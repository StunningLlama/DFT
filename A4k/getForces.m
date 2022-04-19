function F = getForces(W)
global gbl_S; global gbl_R; global gbl_G2; global gbl_X; global gbl_G; global gbl_Sf; global gbl_X; global gbl_f; global gbl_r; global gbl_Z;
global gbl_n;
S = gbl_S;
R = gbl_R;
G2 = gbl_G2;
X = gbl_X;
G = gbl_G;
Sf = gbl_Sf;
X = gbl_X;
r = gbl_r;
Z = gbl_Z;
n = gbl_n;

dr= sqrt(sum((r - ones(prod(S), 1)*sum(R,2)'/2).^2, 2));
sigma1=0.25;
g1=exp(-dr.^2/(2*sigma1^2))/sqrt(2*pi*sigma1^2)^3;

for alpha = 1:size(X, 1)
    for j = 1:3
        dSf = -i*G(:,j).*(exp(-i*G*X(alpha,:)'));
        dVtilde = (-4*pi*Z./(G2)).*dSf;
        dVtilde(1)=0.;
        
        nnuc=cI(cJ(g1).*Sf); nnuc=real(nnuc);
        phi=cI(Linv(-4*pi*O(cJ(nnuc))));
        phi=real(phi);
        
        dn=cI(cJ(g1).*dSf); dn=real(dn);
        dphi=cI(Linv(-4*pi*O(cJ(dn))));
        dphi=real(dphi);
        
        dU=0.5*real(cJ(dphi)'*O(cJ(nnuc))) + 0.5*real(cJ(phi)'*O(cJ(dn)));
        dE(alpha, j) = real(cJ(dVtilde).' * n) + dU;
    end
end
F = dE;
end
