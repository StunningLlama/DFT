function main
    h = 0.001;
    fplot(@(x) (logistic(x+h, 3.4)-logistic(x, 3.4))/h, [0 1]);
end

function y = logistic(x, a)
    y = 0;
    for i = [1:1:20]
    x = a*x.*(1-x);
    y = y+x/factorial(i)*(-1)^(i+1);
    end
end