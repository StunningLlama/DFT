function sums = complexinnerprod(Wa,Wb)
sums = 0;
for k = [1:size(Wa,2)]
    sums = sums + real(sum(sum(conj(Wa{k}).*Wb{k})));
end
end