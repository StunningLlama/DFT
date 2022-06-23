function [out]=pccg(W,Nit, cgform)
W = orthonormalize(W);
alphat = 3e-5;
eprev = getE(W);
for it = 1:1:Nit
    g = getgrad(W);
    
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
    
    gt = getgrad(linadd(W,d,1,alphat));
    alpha = alphat*(complexinnerprod(g,d))/complexinnerprod(linadd(g,gt,1,-1),d);
    W = linadd(W,d,1,alpha);
    gprev = g;
    dprev = d;
    E = getE(W);
    disp(E); %# <= New statements
    gtnorm = getnorm(gt);
    %disp2(gtnorm);
    if (it > 2)
        if (abs(E-eprev) < 1e-14)
            break;
        end
    end
    eprev = E;
%    disp("Ε: " + num2str(Elist(it)));
end
out = W;
end