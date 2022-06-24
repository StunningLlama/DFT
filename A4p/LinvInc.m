function out=LinvInc(in, k, factor)
global gbl_R;
global gbl_G;
global gbl_kvectors;


kvec = gbl_kvectors(k,:);
karray = ones(size(in, 1), 1)*kvec;

G2= sum((gbl_G+factor*karray).^2, 2);


out= -(in./G2)/det(gbl_R); %# <=== YOUR CODE HERE
end