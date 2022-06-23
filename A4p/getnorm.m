function Wnorm = getnorm(W)
normsq = 0;
elements = 0;
for k = [1:size(W,2)]
    normsq = normsq + sum(sum(abs(W{k}).^2));
    elements = elements + prod(size(W{k}));
end
Wnorm = sqrt(normsq)/elements;
end