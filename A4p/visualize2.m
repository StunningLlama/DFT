function visualize2(n)

global gbl_S; global gbl_M; global gbl_Ns; global gbl_kpoints;
S = gbl_S;
M = gbl_M;
Ns = gbl_Ns;
global gbl_R;

%levelsets = [0.1:0.1:0.9];
%alpha = 0.2;
levelsets = [0.1 0.3 0.5 0.7 0.9];
alpha = 0.2;
interval = [0 1 0 1 0 1];
%Psi = W;
dat=n;
bigarray = zeros([S(2) S(1) S(3)]);
bigarray(sub2ind(size(bigarray), M(:,2)+1, M(:,1)+1,M(:,3)+1)) = dat;
maxamplitude = max(bigarray, [], 'all')

for level = levelsets
    fimplicit3(@(xq,yq,zq) (interp3(((S(1):-1:1) - 0.5)/S(1),((S(2):-1:1) - 0.5)/S(2),((S(3):-1:1) - 0.5)/S(3), bigarray-level*maxamplitude, xq, yq, zq)),'EdgeColor','none','FaceAlpha',alpha);
    hold on;
end
hold off;
%    isosurface((1:S(1))/S(1),(1:S(2))/S(2),(1:S(3))/S(3), bigarray, 0.01);
axis(interval);
pause();
end