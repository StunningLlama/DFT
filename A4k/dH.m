function out=dH(Y, Vsp, k)
global gbl_Vdual;
global gbl_f;
%Y=W*inv(sqrtm(W'*O(W)));
%uinv = inv(W'*O(W));
%dUinv= -uinv*(dW'*O(W) + W'*O(dW))*uinv;
%n2 = gbl_f*diagouter(cI(W*uinv), cI(W));
%disp("test")
%disp2(n)
%disp2(n2)

%dn = gbl_f*diagouter(cI(dW*uinv), cI(W))+gbl_f*diagouter(cI(W*dUinv), cI(W))+gbl_f*diagouter(cI(W*uinv), cI(dW));
%disp2(dn)

global gbl_n;
global gbl_excpn;
global gbl_excppn;
global gbl_JdagOJn;
global gbl_IW;
global gbl_IY;

n = gbl_n;
IW = gbl_IW;
out = zeros(size(Y));
for col=1:size(Y,2)
out(:,col) = out(:,col) + cIdag(Diagprod(Vsp, IW(:,col,k)));
end
end