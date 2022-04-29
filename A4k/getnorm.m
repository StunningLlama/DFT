function Wnorm = getnorm(W)
global gbl_kpoints;
global gbl_weights;
normsq = 0;
elements = 0;
for k = [1:gbl_kpoints]
    normsq = normsq + sum(sum(abs(W{k}).^2));
    elements = elements + prod(size(W{k}));
end
Wnorm = sqrt(normsq)/elements;
end