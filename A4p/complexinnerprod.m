function sums = complexinnerprod(Wa,Wb)
global gbl_kpoints;
sums = 0;
for k = [1:gbl_kpoints]
    sums = sums + real(sum(sum(conj(Wa{k}).*Wb{k})));
end
end