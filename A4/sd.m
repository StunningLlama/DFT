function [out, Elist]=sd(W, Nit)
W = W*sqrtm(inv(W'*O(W)));
Elist = zeros(Nit, 1);
out = W;
alpha = 1e-3;
for it = 1:1:Nit
    out = out - alpha*getgrad(out);
    Elist(it)=getE(out); %# <= New statements
    Elist(it)
end
end