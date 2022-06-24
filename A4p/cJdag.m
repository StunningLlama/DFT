function out=cJdag(in)
global gbl_S;
global gbl_gperm;
out = zeros(prod(gbl_S), size(in, 2));
for col = 1:size(in, 2)
    full=zeros(prod(gbl_S),1); full(gbl_gperm)=in(:,col);
    out(:,col) = fft3(full,gbl_S,1)/prod(gbl_S);
end
end