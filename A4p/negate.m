function out=negate(W)
out = {};
for k = [1:size(W,2)]
    out{k} = -W{k};
end
end