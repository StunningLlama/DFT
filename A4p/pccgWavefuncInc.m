% out=pccgWavefunc(W, dTau, Nit, cgform)
function out=pccgWavefuncInc(W, dVext, q, Nit, cgform)

alphat = 3e-5;

errs = [];

b = negate(getPsiVspDerivWFillingsInc(dVext, q));
dW = b;
gtnormprev = 1e99;

for it = 1:1:Nit
    g = linadd(getPsiPsiDerivWFillingsInc(dW, q), b, 1,-1);
    
    if (it > 1)
        disp("Angle cosine: " + num2str(complexinnerprod(g,d)/sqrt(complexinnerprod(g,g)*complexinnerprod(d,d))));
        disp("CG test: " + num2str(complexinnerprod(g,KInc(gprev,q))/sqrt(complexinnerprod(g,KInc(g,q))*complexinnerprod(gprev,KInc(gprev,q)))));
        
        beta = 0;
        if (cgform == 1)
            beta = complexinnerprod(g,KInc(g,q))/complexinnerprod(gprev,KInc(gprev,q));
        elseif (cgform == 2)
            %beta = (real(sumall(conj(g-gprev).*K(g)))/real(sumall(conj(gprev).*K(gprev))));
            disp("NOT IMPLEMENTED");
        elseif (cgform == 3)
            %beta = (real(sumall(conj(g-gprev).*K(g)))/real(sumall(conj(g-gprev).*dprev)));
            disp("NOT IMPLEMENTED");
        end
        d = linadd(KInc(g, q), dprev, -1, beta);
    else
        d = negate(KInc(g, q));
    end
    
    gt = linadd(getPsiPsiDerivWFillingsInc(linadd(dW,d,1,alphat), q), b,1,-1);
    alpha = alphat*(complexinnerprod(g,d))/complexinnerprod(linadd(g,gt,1,-1),d);
    dW = linadd(dW,d,1,alpha);
    gprev = g;
    dprev = d;
    gtnorm = getnorm(mult(d,alpha));
    errs = [errs gtnorm];
    disp2(gtnorm);
    
    if (it > 5 && gtnorm > gtnormprev)
       break;
    end
    gtnormprev = gtnorm;
    out = dW;
end

plot(log(errs));