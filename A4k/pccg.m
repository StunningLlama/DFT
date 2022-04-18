function [out]=pccg(W,Nit, cgform)
W = orthonormalize(W);
alphat = 3e-5;
for it = 1:1:Nit
    g = getgrad(W);
    
    if (it > 1)
%        disp("Angle cosine: " + num2str(real(trace(g'*d))/sqrt(real(trace(g'*g))*real(trace(d'*d)))));
%        disp("CG test: " + num2str(real(trace(g'*K(gprev)))/sqrt(real(trace(g'*K(g)))*real(trace(gprev'*K(gprev))))));
        
        beta = 0;
        if (cgform == 1)
            beta = (real(sumall(conj(g).*K(g)))/real(sumall(conj(gprev).*K(gprev))));
        elseif (cgform == 2)
            beta = (real(sumall(conj(g-gprev).*K(g)))/real(sumall(conj(gprev).*K(gprev))));
        elseif (cgform == 3)
            beta = (real(sumall(conj(g-gprev).*K(g)))/real(sumall(conj(g-gprev).*dprev)));
        end
        d = -K(g) + beta*dprev;
    else
        d = -K(g);
    end
    
    gt = getgrad(W+alphat*d);
    alpha = alphat*(real(sumall(conj(g).*d)))/(real(sumall(conj(g-gt).*d)));
    W = W + alpha*d;
    gprev = g;
    dprev = d;
    disp(getE(W)); %# <= New statements
    gtnorm = sqrt(sumall(abs(gt).^2))/prod(size(gt));
    %disp2(gtnorm);
    if (it > 2)
        if (gtnorm < 1e-10)
            break;
        end
    end
%    disp("Ε: " + num2str(Elist(it)));
end
out = W;
end