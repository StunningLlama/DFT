function out=L(in, k)

global gbl_R
global gbl_Gc
global gbl_kvectors
global gbl_G2c;

if (isempty(k))
    out= -det(gbl_R)*(gbl_G2c{0}*ones(1,size(in,2))).*in;
else
    kvec = gbl_kvectors(k,:);
    karray = ones(size(in, 1), 1)*kvec;
    
    out= -det(gbl_R)*(sum((gbl_Gc{k}+karray).^2, 2)*ones(1,size(in,2))).*in;
end
end