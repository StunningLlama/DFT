function out=sd(W, Nit)
out = W
alpha = 3e-5;
for i = 1:1:Nit
    out = out - alpha*getgrad(out);
    disp(getE(out));
end
end