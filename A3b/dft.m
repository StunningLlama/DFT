global gbl_Vdual;
%# Set the orbital occupancies
global gbl_f;
gbl_f=2; %# The usual case of a constant filling of two electrons per orbital

omega = 2;
V = 0.5*omega^2*sum((r - ones(prod(S), 1)*sum(R,2)'/2).^2, 2);
gbl_Vdual = cJdag(O(cJ(V)));

%# Finite difference test
Ns=4; %# Number of states
randn('seed',0.2004);
W=(randn(prod(S),Ns)+i*randn(prod(S),Ns));
more off; %# View output as it is computed
fdtest(W);
%pause();

%
%W = sd(W, 250);

%# Allow for more digits in printouts
format long
%# Converge using steepest descents

%# Converge
W=sd(W,20); %# 20 iterations of simple sd() to get nearer to the minimum
W=W*inv(sqrtm(W'*O(W))); %# Restart as orthonormal functions
[Wlm, Elm]=lm(W,50); %# 50 iterations of lm() from W, while recording results
[Wpclm, Epclm]=pclm(W,50); %# 50 iterations of pclm from same W
[Wcg1, Ecg1]=pccg(W,50,1); %# 50 iters from W of cg versions 1, 2, 3 ...
[Wcg2, Ecg2]=pccg(W,50,2);
[Wcg3, Ecg3]=pccg(W,50,3);
%# Plot results on log scale with nice labels
grid on;
xlabel("Iteration (#) ->"); ylabel("Error (H)");
semilogy([1:length(Elm)],Elm-43.3371147782040,'-', ...
[1:length(Epclm)],Epclm-43.3371147782040,'-', ...
[1:length(Ecg1)],Ecg1-43.3371147782040,'-', ...
[1:length(Ecg2)],Ecg2-43.3371147782040,'-', ...
[1:length(Ecg3)],Ecg3-43.3371147782040,'-');
legend('lm', 'pclm', 'pccg-FR', 'pccg-PR', 'pccg-HS');
%# Set W to best result for plotting
W=Wcg3;
pause; %# Wait for keyboard input before continuing


%# Extract and display final results
[Psi, epsilon]=getPsi(W);
for st=1:Ns
fprintf('=== State # %d, Energy = %f ===\n',st,epsilon(st));
dat=abs(cI(Psi(:,st))).^2;
for k=1:3
sl=slice(dat,S,S(k)/2,k);
name=sprintf("psi%d_m%d.ppm",st,k);
%ppm(name,sl*0.3,sl,sl);% system(["display " name "&"]);
contourf(sl);
pause();
end
end