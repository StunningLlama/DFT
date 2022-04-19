test = 1
if (test==0) % Test visualization
    X=[8 8 8; 8+1.5 8 8];
    setup(X, 2, [2 2], true);
    [W,E1] = iterate(20);
    visualize(W);
elseif (test == 1) % Test atomic force calculation
    X=[8 8 8; 8+2 8 8];
    setup(X, 1, 1, true);
    W1 = iterate(20);
    setupPccgWavefunc(W1);
    E1 = getE(W1) + ewald();
    F = getForces(W1);
    ds=0.01;
    X=[8 8 8; 8+2+ds 8 8]
    setup(X, 1, 1, true);
    W2 = iterate(20);
    E2 = getE(W2) + ewald();
    F
    (E2-E1)/ds
elseif (test == 2) % Direct FD test of dPsiPsi
    X=[8 8 8; 8+2 8 8];
    setup(X, 1, 1, true);
    %[W,E1] = iterate(20);
    rng('default');
    rng(245);
    global gbl_active;
    global gbl_Ns;
    global gbl_kpoints;
    W=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dWa=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dWb=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dWa = dWa*0.00001;
    dWb = dWb*0.00001;
    
    %dWa = dWa*0.000000;
    %dWb = dWb*0.000000;
    %dWa(3,1)=0.0001;
    %dWb(2,1)=0.0001;
    
    deltaE = (getE(W+dWa+dWb) - getE(W+dWa)) - (getE(W+dWb) - getE(W));
    deltaE
    
    setupPccgWavefunc(W);
    dwGradE = getPsiPsiDerivWFillings(dWb);
    2*real(sumall(conj(dWa).*dwGradE))
    %dGrad = getgrad(W+dWa) - getgrad(W);
    %dGrad2 = getPsiPsiDeriv(W, dWa);
    
    %A = H(W+dWa)-H(W);
    %B = dH(W, dWa) + H2(W, dWa);
    %Bp = H(dWa);
elseif (test == 2.5) % Direct FD test of dPsiPsi
    X=[8 8 8; 8+2 8 8];
    setup(X, 1, 1, true);
    %[W,E1] = iterate(20);
    rng('default');
    rng(245);
    global gbl_active;
    global gbl_Ns;
    global gbl_kpoints;
    W=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dWa=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dWa = dWa*0.00000001;
    %dWb = dWb*0.0000001;
    
    %dWa = dWa*0.000000;
    %dWb = dWb*0.000000;
    %dWa(3,1)=0.0001;
    %dWb(2,1)=0.0001;
    
    %dwGradE = getPsiPsiDeriv(W, dWb);
    %2*real(trace(dWa'*dwGradE))
    dGrad = getgrad(W+dWa) - getgrad(W);
    setupPccgWavefunc(W);
    disp("Phase 2:");
    dGrad2 = getPsiPsiDerivWFillings(dWa);
    %dGrad3 = getPsiPsiDeriv(W, dWa);
    disp2(dGrad);
    disp2(dGrad2);
    %disp2(dGrad3);
    %A = H(W+dWa)-H(W);
    %B = dH(W, dWa) + H2(W, dWa);
    %Bp = H(dWa);
elseif (test==3) %FD test dTau
    global gbl_active;
    global gbl_Ns;
    global gbl_kpoints;
    X=[8 8 8; 8+2 8 8];
    rng('default');
    rng(254);
    dX=(randn(size(X)));
    dX = dX*0.0000001;
    
    setup(X, 1, 1, true);
    W=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    setupPccgWavefunc(W);
    Grad1 = getgrad(W);
    setup(X+dX, 1, 1, true);
    Grad2 = getgrad(W);
    dGrad = Grad2-Grad1;
    
    setup(X, 1, 1, true);
    dGradAnal = getPsiTauDerivWFillings(dX);
    %dGradAnal2 = getPsiTauDeriv(W, dX);
elseif (test==4)
    X=[8 8 8; 8+2 8 8];
    setup(X, 1, 1, true);
    %[W,E1] = iterate(20);
    rng('default');
    rng(249);
    global gbl_active;
    global gbl_Ns;
    global gbl_kpoints;
    W=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dW=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dW = dW*0.01;
    dWa = dW*0;
    dWb = dW*0;
    
    indexA = 380;
    indexB = 1020;
    dWa(indexA) = 1;
    dWb(indexB) = 1;
    
    dGrad1 = getPsiPsiDeriv(W, dW);
    dGrad1 = getPsiPsiDeriv(W, dWa);
    c1 = dGrad1(indexB);
    dGrad2 = getPsiPsiDeriv(W, dWb);
    c2 = dGrad2(indexA);
    
    c1
    c2
    disp2(dW'*dGrad1)
    %dGrad = getgrad(W+dWa) - getgrad(W);
    %dGrad2 = getPsiPsiDeriv(W, dWa);
    
    %A = H(W+dWa)-H(W);
    %B = dH(W, dWa) + H2(W, dWa);
    %Bp = H(dWa);
elseif (test==5) %Test that dPsiPsi is symmetric
    
    X=[8 8 8; 8+2 8 8];
    setup(X, 1, 1, true);
    %[W,E1] = iterate(20);
    rng('default');
    rng(249);
    global gbl_active;
    global gbl_Ns;
    global gbl_kpoints;
    W=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    dW=(randn(2*length(gbl_active),gbl_Ns, gbl_kpoints));
    %dW = dW*0.01;
    dWa = dW*0;
    dWb = dW*0;
    
    indexA = 3800;
    indexB = 29000;
    dWa(indexA) = 1;
    dWb(indexB) = 1;
    
    %dGrad1 = getdGradtmp(W, dW);
    dGrad1 = getdGradtmp(W, dWa);
    c1 = dGrad1(indexB);
    dGrad2 = getdGradtmp(W, dWb);
    c2 = dGrad2(indexA);
    
    c1
    c2
    %disp2(dW'*dGrad1)
elseif (test==6) % Test conjugate gradient for perturbation
    X=[8 8 8; 8+2 8 8];
    setup(X, 1, 1, true);
    
    dX=(randn(size(X)));
    
    [W,E1] = iterate(20);
    disp("Done pt. 1");
    setupPccgWavefunc(W);
    dW = 0.001*pccgWavefunc(W, dX, 100, 1);
    
    
    grad1 = getgrad(W);
    setup(X+dX*0.001, 1, 1, true);
    grad2 = getgrad(W+dW);
    disp2(grad1);
    disp2(grad2);

elseif (test == 7) % Direct FD test of dTauTau
    X=[8 8 8; 8+2 8 8];
    rng('default');
    rng(259);
    dXa=(randn(size(X)));
    dXa = dXa*0.0001;
    dXb=(randn(size(X)));
    dXb = dXb*0.0001;
    
    setup(X, 1, 1, true);
    %[W,E1] = iterate(20);
    global gbl_active;
    global gbl_Ns;
    global gbl_kpoints;
    W=(randn(length(gbl_active),gbl_Ns, gbl_kpoints)+i*randn(length(gbl_active),gbl_Ns, gbl_kpoints));
    
    setupPccgWavefunc(W);
    dxGradE = getTauTauDeriv(dXa, dXb)
    
    E1=getE(W);
    setup(X+dXa, 1, 1, true);
    E2=getE(W);
    setup(X+dXb, 1, 1, true);
    E3=getE(W);
    setup(X+dXa+dXb, 1, 1, true);
    E4=getE(W);
    deltaE = E4-E3-E2+E1;
    deltaE
    
    %dWa = dWa*0.000000;
    %dWb = dWb*0.000000;
    %dWa(3,1)=0.0001;
    %dWb(2,1)=0.0001;
    
    %dGrad = getgrad(W+dWa) - getgrad(W);
    %dGrad2 = getPsiPsiDeriv(W, dWa);
    
    %A = H(W+dWa)-H(W);
    %B = dH(W, dWa) + H2(W, dWa);
    %Bp = H(dWa);
elseif(test == 8)
    X=[8 8 8; 8+2 8 8];
    rng('default');
    rng(259);
    dXa=zeros(size(X));
    dXb=zeros(size(X));
    dXa(1, 1) = 1;
    dXb(2, 1) = 1;
    setup(X, 1, 1, true);
    [W,E1] = iterate(20);
    setupPccgWavefunc(W);
    dWa = pccgWavefunc(W, dXa, 100, 1);
    dWb = pccgWavefunc(W, dXb, 100, 1);
    Kanal = calcSpringConstant(W, X, dXa, dXb, dWa, dWb)
    h = 0.001;
    
    setup(X+h*dXa, 1, 1, true);
    W2=pccg(W,50,1);
    E2=getE(W2);
    
    setup(X+h*dXb, 1, 1, true);
    W3=pccg(W,50,1);
    E3=getE(W3);
    
    setup(X+h*dXa+h*dXb, 1, 1, true);
    W4=pccg(W,50,1);
    E4=getE(W4);
    
    Knum = -(E4-E3-E2+E1)/h^2
elseif (test==9)
    Es = []
    for dx = [-0.5:0.1:0.5]
        X = [0 0 0; 2+dx 0 0];
        setup(X, 1, 1, true);
        [W,E1] = iterate(20);
        Es = [Es getE(W)];
    end
    Es
elseif (test == 10) % Test of structure factor 2nd deriv
    X=[8 8 8; 8+2 8 8];
    rng('default');
    rng(254);
    dXa=(randn(size(X)));
    dXa = dXa*0.0001;
    dXb=(randn(size(X)));
    dXb = dXb*0.0001;
    
    setup(X, 1, 1, true);
    %[W,E1] = iterate(20);
    Sf1=getSf(X);
    Sf2=getSf(X+dXa);
    Sf3=getSf(X+dXb);
    Sf4=getSf(X+dXa+dXb);
    deltaSf = Sf4-Sf3-Sf2+Sf1;
    disp2(deltaSf);
    dSfAnal = getdsqSf(X, dXa, dXb);
    disp2(dSfAnal);
    
  
    
    setup(X, 1, 1, true);
    E1=ewald();
    setup(X+dXa, 1, 1, true);
    E2=ewald();
    setup(X+dXb, 1, 1, true);
    E3=ewald();
    setup(X+dXa+dXb, 1, 1, true);
    E4=ewald();
    deltaE = E4-E3-E2+E1
    dEAnal = getdsqEwald(X, dXa, dXb)
    
    %dWa = dWa*0.000000;
    %dWb = dWb*0.000000;
    %dWa(3,1)=0.0001;
    %dWb(2,1)=0.0001;
    
    %dGrad = getgrad(W+dWa) - getgrad(W);
    %dGrad2 = getPsiPsiDeriv(W, dWa);
    
    %A = H(W+dWa)-H(W);
    %B = dH(W, dWa) + H2(W, dWa);
    %Bp = H(dWa);
elseif (test == 11) %Combined spring const test w/ electronic and nuclear energy
    
    X=[8 8 8; 8+1.45 8 8];
    rng('default');
    rng(259);
    dXa=zeros(size(X));
    dXb=zeros(size(X));
    dXa(1, 1) = 1;
    dXb(2, 1) = 1;
    setup(X, 1, 1, true);
    W = iterate(20);
    E1 = getE(W)+ewald();
    
    setupPccgWavefunc(W);
    dWa = pccgWavefunc(W, dXa, 100, 1);
    dWb = pccgWavefunc(W, dXb, 100, 1);
    Kanal = calcSpringConstant(W, X, dXa, dXb, dWa, dWb)-getdsqEwald(X, dXa, dXb);
    h = 0.001;
    
    setup(X+h*dXa, 1, 1, true);
    W2=pccg(W,50,1);
    E2=getE(W2)+ewald();
    
    setup(X+h*dXb, 1, 1, true);
    W3=pccg(W,50,1);
    E3=getE(W3)+ewald();
    
    setup(X+h*dXa+h*dXb, 1, 1, true);
    W4=pccg(W,50,1);
    E4=getE(W4)+ewald();
    
    Knum = -(E4-E3-E2+E1)/h^2
end
%visualize(W);