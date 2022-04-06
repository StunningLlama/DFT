function setupPccgWavefunc(W)
global gbl_f;
global gbl_U;
global gbl_Uinv;
global gbl_Usqrtinv;
global gbl_Htilde;
global gbl_Y;
global gbl_tmp1;
global gbl_UsiFUsi;
global gbl_WUsiFUsi;
global gbl_HFmFH;
global gbl_QHFmFH;
global gbl_HWUsi;
global gbl_n;
global gbl_excpn;
global gbl_excppn;
global gbl_hsum;
global gbl_Vdual;
global gbl_JdagOJn;
global gbl_IW;
global gbl_IY;

gbl_U = W'*O(W);
gbl_Uinv = inv(gbl_U);
gbl_Usqrtinv = sqrtm(gbl_Uinv);
gbl_Htilde = gbl_Usqrtinv*(W'*H(W))*gbl_Usqrtinv;
gbl_Y=W*gbl_Usqrtinv;
gbl_tmp1 = H(W*gbl_Usqrtinv*diag(gbl_f)*gbl_Usqrtinv);
gbl_UsiFUsi = gbl_Usqrtinv*diag(gbl_f)*gbl_Usqrtinv;
gbl_WUsiFUsi = W*gbl_UsiFUsi;
gbl_HFmFH = gbl_Htilde*diag(gbl_f)-diag(gbl_f)*gbl_Htilde;
gbl_QHFmFH = Q(gbl_Htilde*diag(gbl_f)-diag(gbl_f)*gbl_Htilde, gbl_U);
gbl_HWUsi = H(W)*gbl_Usqrtinv;

gbl_n = getn(cI(gbl_Y), gbl_f);
gbl_excpn = excpVWN(gbl_n);
gbl_excppn = excppVWN(gbl_n);
gbl_hsum = gbl_Vdual + cJdag(O(-4*pi*Linv(O(cJ(gbl_n))))) ...
    + cJdag(O(cJ(excVWN(gbl_n)))) ...
    + Diagprod(excpVWN(gbl_n), cJdag(O(cJ(gbl_n))));
gbl_JdagOJn = cJdag(O(cJ(gbl_n)));

gbl_IW = cI(W);
gbl_IY = cI(gbl_Y);
end