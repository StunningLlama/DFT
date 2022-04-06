function out=cIdag(in)
global gbl_S;
global gbl_active;
out=zeros(length(gbl_active),size(in,2));
for col = 1:size(in, 2)
    full=fft3(in(:,col),gbl_S,-1); out(:,col)=full(gbl_active);
end
end