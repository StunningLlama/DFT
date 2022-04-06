function out=K(W)
global gbl_G2;
global gbl_G2c
if size(W,1)==length(gbl_G2c)
    out= (1./(1+gbl_G2c)*ones(1,size(W,2))).*W;
else
    out= (1./(1+gbl_G2)*ones(1,size(W,2))).*W;
end
end