function out=cI(in)
global gbl_S;
out = fft3(in,gbl_S,1);
end