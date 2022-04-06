function out=K(W)
global gbl_G2;
out= (1./(1+gbl_G2)*ones(1,size(W,2))).*W;
end