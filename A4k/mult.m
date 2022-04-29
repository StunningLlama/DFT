function out=mult(W,a)
global gbl_kpoints;
out = {};
for k = [1:gbl_kpoints]
    out{k} = W{k}*a;
end
end