% out=pccgWavefunc(W, dTau, Nit, cgform)
function out=pccgWavefunc(W, dTau, Nit, cgform)

%W = W*sqrtm(inv(W'*O(W)));
alphat = 3e-5;

b = -getPsiTauDerivWFillings(W, dTau);
dW = b;
for it = 1:1:Nit
    g = getPsiPsiDerivWFillings(W, dW) - b;
    
    if (it > 1)
%        disp("Angle cosine: " + num2str(real(trace(g'*d))/sqrt(real(trace(g'*g))*real(trace(d'*d)))));
%        disp("CG test: " + num2str(real(trace(g'*K(gprev)))/sqrt(real(trace(g'*K(g)))*real(trace(gprev'*K(gprev))))));
        
        beta = 0;
        if (cgform == 1)
            beta = (real(trace(g'*K(g)))/real(trace(gprev'*K(gprev))));
        elseif (cgform == 2)
            beta = (real(trace((g-gprev)'*K(g)))/real(trace(gprev'*K(gprev))));
        elseif (cgform == 3)
            beta = (real(trace((g-gprev)'*K(g)))/real(trace((g-gprev)'*dprev)));
        end
        d = -K(g) + beta*dprev;
    else
        d = -K(g);
    end
    
    gt = getPsiPsiDerivWFillings(W, dW+alphat*d) - b;
    alpha = alphat*(real(trace(g'*d)))/(real(trace((g-gt)'*d)));
    dW = dW + alpha*d;
    gprev = g;
    dprev = d;
    gtnorm = norm(gt)/prod(size(gt));
    disp2(gtnorm);
    if (it > 10)
        if (gtnorm < 1e-7)
            break;
        end
    end
%    disp("Ε: " + num2str(Elist(it)));
out = dW;
end