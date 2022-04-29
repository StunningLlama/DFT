function out=sd(W, Nit)
W = orthonormalize(W);
out = W;
alpha = 0.001;
%Eprev = 10000;
for it = 1:1:Nit
    alphahigh = alpha*1.2;
    alphalow = alpha/1.2;
    outhigh = linadd(out,K(getgrad(out)),1,-alphahigh);
    outlow = linadd(out,K(getgrad(out)),1,-alphalow);
    Ehigh = getE(outhigh);
    Elow = getE(outlow);
    if (Elow < Ehigh)
        out = outlow;
        alpha = alphalow;
        Elow
    else
        out = outhigh;
        alpha = alphahigh;
        Ehigh
    end
    alpha
    
    
    
%     out = out - alpha*K(getgrad(out));
%     E = getE(K(out))
    %disp("Doing something");
    %Eprev = E;
end
out = K(out);
end