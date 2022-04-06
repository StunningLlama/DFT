function out=cI(in)
global gbl_active;
global gbl_S;
out = zeros(prod(gbl_S), size(in, 2));
for col = 1:size(in, 2)
    if size(in,1)==prod(gbl_S)
        out(:,col) = fft3(in(:,col),gbl_S,1);
    else
        full=zeros(prod(gbl_S),1); full(gbl_active)=in(:,col);
        out(:,col)= fft3(full,gbl_S,1); %# previous code with "full" replacing "in(:,col)"
    end
end
end