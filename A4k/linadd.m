function out=linadd(Wa,Wb,a,b)
global gbl_kpoints;
out = {};
for k = [1:gbl_kpoints]
    out{k} = Wa{k}*a+Wb{k}*b;
end
end