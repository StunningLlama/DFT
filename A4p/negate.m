function out=negate(W)
global gbl_kpoints;
out = {};
for k = [1:gbl_kpoints]
    out{k} = -W{k};
end
end