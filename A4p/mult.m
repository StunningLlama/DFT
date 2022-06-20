function out=mult(W,a)
out = {};
for k = [1:size(W,2)]
    out{k} = W{k}*a;
end
end