function setupPccgWavefuncInc(W)
global gbl_f;
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
global gbl_kpoints;
global gbl_weights;
global gbl_W;

gbl_Htilde = {};
gbl_Y = {};
gbl_tmp1 = {};
gbl_UsiFUsi = {};
gbl_WUsiFUsi = {};
gbl_HFmFH = {};
gbl_QHFmFH = {};
gbl_HWUsi = {};

gbl_W = W;
gbl_Y = W;

gbl_IW = {};
gbl_IY = {};
for k = [1:gbl_kpoints]
    gbl_IW{k} = cI(W{k}, k);
    gbl_IY{k} = gbl_IW{k};
end

gbl_n = getn(gbl_Y, gbl_f);

gbl_excpn = excpVWN(gbl_n);
gbl_excppn = excppVWN(gbl_n);
gbl_hsum = gbl_Vdual + cJdag(O(-4*pi*Linv(O(cJ(gbl_n))))) ...
    + cJdag(O(cJ(excVWN(gbl_n)))) ...
    + Diagprod(excpVWN(gbl_n), cJdag(O(cJ(gbl_n))));
gbl_JdagOJn = cJdag(O(cJ(gbl_n)));

for k = [1:gbl_kpoints]
    gbl_Htilde{k} = W{k}'*H(W{k}, gbl_hsum, k);
    gbl_tmp1{k} = H(W{k}*diag(gbl_f), gbl_hsum, k);
    gbl_UsiFUsi{k} = diag(gbl_f);
    gbl_WUsiFUsi{k} = W{k}*gbl_UsiFUsi{k};
    gbl_HFmFH{k} = gbl_Htilde{k}*diag(gbl_f)-diag(gbl_f)*gbl_Htilde{k};
    gbl_QHFmFH{k} = (gbl_Htilde{k}*diag(gbl_f)-diag(gbl_f)*gbl_Htilde{k})/2.0;
    gbl_HWUsi{k} = H(W{k}, gbl_hsum, k);
end
end