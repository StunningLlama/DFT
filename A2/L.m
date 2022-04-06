
function out=L(in)
global gbl_R;
global gbl_G2;
out= -det(gbl_R)*(gbl_G2*ones(1,size(in,2))).*in; %# <=== YOUR CODE HERE
end