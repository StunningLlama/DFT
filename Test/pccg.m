function [out, Elist]=pccg(W,Nit, cgform)
%W = W*sqrtm(inv(W'*O(W)));
Elist = zeros(Nit, 1);
alphat = 3e-5;
for it = 1:1:Nit
    g = getgrad(W);
    
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
    
    gt = getgrad(W+alphat*d);
    alpha = alphat*(real(trace(g'*d)))/(real(trace((g-gt)'*d)));
    W = W + alpha*d;
    gprev = g;
    dprev = d;
    Elist(it)=getE(W); %# <= New statements
    Elist(it)
    if (it > 10)
        if (abs((Elist(it) - Elist(it-1))/Elist(it)) < 0.0000000001)
            Elist(Nit) = Elist(it);
            break;
        end
    end
%    disp("Ε: " + num2str(Elist(it)));
end
out = W;
end