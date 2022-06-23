% Calculates differential of grad E given differential dTau.
function dwGradE = getPsiTauDerivWFillings(dTau)

global gbl_f;
global gbl_Y;
global gbl_kpoints;
global gbl_weights;

F = diag(gbl_f);

Y = gbl_Y;

dHtilde = {};
dwGradE = {};

for k = [1:gbl_kpoints]
    dHtilde{k} = Y{k}'*dHtau(Y{k}, dTau, k);
    
    tmp = dHtau(Y{k}*F, dTau, k);
    dwGradE{k} = tmp-O(Y{k}*(Y{k}'*tmp)) + O(Y{k}*(dHtilde{k}*F-F*dHtilde{k})/2.0);
    dwGradE{k} = dwGradE{k}*gbl_weights(k);
end
end