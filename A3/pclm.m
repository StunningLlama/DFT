function [out, Elist]=pclm(W,Nit)
W = W*sqrtm(inv(W'*O(W)));
Elist = zeros(Nit, 1);
alphat = 3e-5;
for it = 1:1:Nit
g = getgrad(W);
if (it > 1)
    disp("Angle cosine: " + num2str(real(trace(g'*d))/sqrt(real(trace(g'*g))*real(trace(d'*d)))))
end
d = -K(g);
gt = getgrad(W+alphat*d);
alpha = alphat*(real(trace(g'*d)))/(real(trace((g-gt)'*d)));
W = W + alpha*d;
Elist(it)=getE(W); %# <= New statements
Elist(it);
disp(Elist(it));
end
out = W;
end