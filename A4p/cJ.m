function out=cJ(in)
global gbl_S;
global gbl_gperm;
out=zeros(length(gbl_gperm),size(in,2));
for col = 1:size(in, 2)
    full=fft3(in(:,col),gbl_S,-1)/prod(gbl_S);
    out(:,col)=full(gbl_gperm);
end
end