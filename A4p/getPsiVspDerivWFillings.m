% Calculates differential of grad E given differential dTau.
function dwGradE = getPsiVspDerivWFillings(dVsp)

global gbl_f;
global gbl_Y;
global gbl_kpoints;
global gbl_weights;

F = diag(gbl_f);

Y = gbl_Y;

dHtilde = {};
dwGradE = {};

for k = [1:gbl_kpoints]
    dHtilde{k} = Y{k}'*dH(cI(Y{k},k), dVsp, k);
    
    tmp = dH(cI(Y{k},k), dVsp, k);
    dwGradE{k} = tmp-O(Y{k}*(Y{k}'*tmp)) + O(Y{k}*(dHtilde{k}*F-F*dHtilde{k})/2.0);
    dwGradE{k} = dwGradE{k}*gbl_weights(k);
end
end