function out=L(in)
global gbl_R;
global gbl_G2;
global gbl_G2c
if size(in,1)==length(gbl_G2c)
    out= -det(gbl_R)*(gbl_G2c*ones(1,size(in,2))).*in;
else
    out= -det(gbl_R)*(gbl_G2*ones(1,size(in,2))).*in;
end
end