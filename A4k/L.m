function out=L(in, k)

global gbl_R
global gbl_Gc
global gbl_kvectors

if (isempty(k))
    out= -det(gbl_R)*(gbl_G2c*ones(1,size(in,2))).*in;
else
    kvec = gbl_kvectors(k,:);
    karray = ones(size(in, 1), 1)*kvec;
    
    out= -det(gbl_R)*(sum((gbl_Gc+karray).^2, 2)*ones(1,size(in,2))).*in;
end
end