function hw4
problem1();
end

function problem1
result = zeros(1,6);
names = ["Saddle", "Stab node", "Unstab node", "Stab spiral", "Unstab spiral", "Other"];
for i = [1:1:100000]
    %type = classify(rand(2)*2-1);
    type = classify(randn(2));
    result(type) = result(type)+1;
end
names
result/sum(result)
end

function type=classify(A)
evs = eigs(A);
if (isreal(evs))
    if (evs(1)<0 & evs(2)<0)
        type = 2;
        return;
    elseif (evs(1) > 0 & evs(2) > 0)
        type = 3;
        return
    elseif (evs(1)*evs(2) < 0)
        type = 1;
        return;
    else
        type = 6;
        return;
    end
else
    if (real(evs(1)) < 0)
        type = 4;
        return;
    elseif (real(evs(1))>0)
        type = 5;
        return;
    else
        type = 6;
        return;
    end
end
end