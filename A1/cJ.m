function out=cJ(in)
global gbl_S;
out = fft3(in,gbl_S,-1)/prod(gbl_S);
end