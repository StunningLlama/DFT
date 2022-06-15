function out=dH(IW, dVsp, k)
global gbl_active;

out=zeros(length(gbl_active{k}),size(IW,2))
for col=1:size(IW,2)
out(:,col) = out(:,col) + cIdag(Diagprod(dVsp, IW(:,col)), k);
end
end