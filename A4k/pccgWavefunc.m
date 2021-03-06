% out=pccgWavefunc(W, dTau, Nit, cgform)
function out=pccgWavefunc(W, dTau, Nit, cgform)

%W = W*sqrtm(inv(W'*O(W)));
alphat = 3e-5;

b = negate(getPsiTauDerivWFillings(dTau));
dW = b;
for it = 1:1:Nit
    g = linadd(getPsiPsiDerivWFillings(dW), b, 1,-1);
    
    if (it > 1)
%        disp("Angle cosine: " + num2str(real(trace(g'*d))/sqrt(real(trace(g'*g))*real(trace(d'*d)))));
%        disp("CG test: " + num2str(real(trace(g'*K(gprev)))/sqrt(real(trace(g'*K(g)))*real(trace(gprev'*K(gprev))))));
        
        beta = 0;
        if (cgform == 1)
            beta = complexinnerprod(g,K(g))/complexinnerprod(gprev,K(gprev));
        elseif (cgform == 2)
            %beta = (real(sumall(conj(g-gprev).*K(g)))/real(sumall(conj(gprev).*K(gprev))));
            disp("NOT IMPLEMENTED");
        elseif (cgform == 3)
            %beta = (real(sumall(conj(g-gprev).*K(g)))/real(sumall(conj(g-gprev).*dprev)));
            disp("NOT IMPLEMENTED");
        end
        d = linadd(K(g), dprev, -1, beta);
    else
        d = negate(K(g));
    end
    
    gt = linadd(getPsiPsiDerivWFillings(linadd(dW,d,1,alphat)), b,1,-1);
    alpha = alphat*(complexinnerprod(g,d))/complexinnerprod(linadd(g,gt,1,-1),d);
    dW = linadd(dW,d,1,alpha);
    gprev = g;
    dprev = d;
    gtnorm = getnorm(mult(d,alpha));
    disp2(gtnorm);
    if (it > 10)
        if (gtnorm < 1e-13)
            break;
        end
    end
%    disp("Ε: " + num2str(Elist(it)));
out = dW;
end