function out=cJdag(in)
global gbl_S;
out = zeros(size(in));
for col = 1:size(in, 2)
    out(:,col) = fft3(in(:,col),gbl_S,1)/prod(gbl_S);
end
end