function out=cJcomp(in, k)
global gbl_S;
global gbl_active;
out = zeros(length(gbl_active{k}),size(in,2));
for col = 1:size(in, 2)
    full=fft3(in(:,col),gbl_S,-1)/prod(gbl_S);
    out(:,col)=full(gbl_active{k});
end
end