function out=linadd(Wa,Wb,a,b)
out = {};
for k = [1:size(Wa,2)]
    out{k} = Wa{k}*a+Wb{k}*b;
end
end