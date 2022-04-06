function main()
%problem1();
%problem3();
problem5();
%problem6();
end

function saveimg(name)
saveas(gcf, "C:\Users\Brandon\Desktop\2218dia\" + name);
end

function problem1()
     x = [-8:0.2:8];
     y1 = exp(x);
     y2 = cos(x);
     figure(1);
     plot(x, y1, x, y2);
     title('Plot of exponential and cosine functions (a)');
     legend('exp(x)', 'cos(x)');
     ylabel('y');
     xlabel('x');
     axis([-Inf, Inf, -5 5]);
     saveimg("4210_1_1.png");
     figure(2);
     for x0 = [-12:1:2]
         [t, x2] = ode45(@(t, x)(exp(x)-cos(x)), [0 5], x0);
         plot(t, x2);
         hold on;
     end
     hold off;
     title("Representative solutions of ODE (d)");
     xlabel("t");
     ylabel("x");
     saveimg("4210_1_2.png");
end

function problem3()
syms k1 a x kinv f fp stab;
f = k1*a*x-kinv*x^2;
roots = solve(f == 0);
fp = diff(f, x);
stab = subs(fp, x, roots);
roots
stab
end

function problem5()
syms x(t) t c
ode = diff(x,t) == -x^c;
initcond = x(0) == 1;
sol = dsolve(ode, initcond);
sol
fplot(subs(sol, c, -0.5));
axis([0 1 0 1]);
title('Plot of solution for c=-1/2');
xlabel('t');
ylabel('x');
saveimg("4210_5_1.png");
end

function problem6()
    cobwebtest('3*x-(x^3)', [-3 3], 1.9, 1, 100, "Cobweb diagram, x0 = 1.9", 1, "xn vs. n, x0 = 1.9");
    cobwebtest('3*x-(x^3)', [-3 3], 2.1, 1, 100, "Cobweb diagram, x0 = 2.1", 2, "xn vs. n, x0 = 2.1");
end

function x = cobwebtest(fcnName, domain, initVal, nStart, nStop, name, num, title2)
clf;

% Get input data (if any).
if (nargin == 0)
    fcnName = '3*x-(x^3)';
    domain  = [-3 3];
    initVal = 2.1;
    nStart  = 1;
    nStop   = 100;
end

% Non-zero only for pedagogical reasons.
snoozeTime = 0.00;

% Start the clock...
tic

% Compute the iterates.
f    = inline(vectorize(fcnName));
x    = zeros(length(initVal),nStop);
x    = initVal;
for n = 1:nStop-1
    x(n+1,:) = f(x(n,:));
end

% Set up the graphics.
h=figure(1);                    % Without these two lines, a MATLAB 
set(h,'doublebuffer','on');     % bug gives an annoying flicker.

xMin = domain(1); xMax = domain(2);
xGrid = linspace(xMin, xMax, 777);  
plot(xGrid,f(xGrid),'r-',...
     'LineWidth',2);            % Plot the user-specified function.
hold on
plot(xGrid, xGrid,'k-');        % Plot the diagonal.
axis([xMin xMax xMin xMax]);    % Scale the picture, and 
axis('square');                 % ensure the plot is square.
title(strcat(['f(x) = ' fcnName '.'], ...
             [' Showing iterates ' num2str(nStart)], ...
             [' to ' num2str(nStop) '.']), 'FontWeight','bold');
wBar = waitbar(0, 'Plotting the iterates...');

for n = nStart:nStop
    waitbar((n - nStart)/(nStop - nStart))
    for k = 1:length(initVal)
        if (n > 1) 
            xplot = [x(n-1,k) x(n,k)];
            yplot = [x(n,k) x(n,k)];
            plot(xplot,yplot,'b')  % Draw horizontal part of cobweb.
            pause(snoozeTime)
        end
        if (n < nStop)
            xplot = [x(n,k) x(n,k)];
            yplot = [x(n,k) x(n+1,k)];
            plot(xplot,yplot,'b')   % Draw vertical part of cobweb.
            pause(snoozeTime)
        end
    end
end 
close(wBar);
disp(['Elapsed time: ' num2str(toc) ' seconds.'])
hold off
title(name);
saveimg("4210_6_" + num2str(num) + "a.png");

figure(2)
plot(nStart:nStop,x,'-d')
title(title2);
xlabel('n');
ylabel('xn');
saveimg("4210_6_" + num2str(num) + "b.png");
end