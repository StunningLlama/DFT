function [out, Elist]=sd(W, Nit)
W = orthonormalize(W);
Elist = zeros(Nit, 1);
out = W;
alpha = 1.5e-3;
for it = 1:1:Nit
    out = out - alpha*getgrad(out);
    disp("Doing something");
    Elist(it)=getE(out); %# <= New statements
    Elist(it)
end
end