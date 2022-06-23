function testpccg()
alphat = 3e-7;

b = [1;1;1;1];

A = diag([1;2;3;4])

x = b;

errs = [];

for it = 1:1:50
    g = A*x-b;
    
    if (it > 1)
        disp("Angle cosine: " + num2str((g'*d)/sqrt((g'*g)*(d'*d))));
        disp("CG test: " + num2str((g'*Ktest(gprev))/sqrt((g'*Ktest(g))*(gprev'*Ktest(gprev)))));
        
        beta = (g'*Ktest(g))/(gprev'*Ktest(gprev));
        d = Ktest(g) - beta*dprev;
    else
        d = -Ktest(g);
    end
    
    gt = A*(x+d*alphat) - b;
    alpha = alphat*(g'*d)/((g-gt)'*d);
    x = x+d*alpha;
    gprev = g;
    dprev = d;
    gtnorm = norm(d*alpha);
    disp2(gtnorm);
    disp2(x);
    pause();
    
    
    errs = [errs gtnorm];
    if (it > 10)
        if (gtnorm < 1e-13)
            break;
        end
    end
    out = x;
end

plot(log(errs));
end

function y = Ktest(x)
%     Ainv = diag([1;1/2;1/3;1/4]);
    Ainv = eye(4);
    y = Ainv*x;
end