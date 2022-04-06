function out=Linv(in)
global gbl_R;
global gbl_G2;

out= -(in./gbl_G2)/det(gbl_R); %# <=== YOUR CODE HERE
out(1) = 0;
end