function visualize(W, X)

global gbl_S; global gbl_M; global gbl_Ns; global gbl_kpoints;
S = gbl_S;
M = gbl_M;
Ns = gbl_Ns;

[Psi, epsilon]=getPsi(W);
for k = [1:gbl_kpoints]
    for st=1:Ns
        dat=abs(cI(Psi(:,st,k))).^2;
        bigarray = zeros([S(1) S(2) S(3)]);
        bigarray(sub2ind(size(bigarray), M(:,1)+1,M(:,2)+1,M(:,3)+1)) = dat;
        maxamplitude = max(bigarray, [], 'all')
        
        %for
        for level = [0.1:0.1:0.9]
            fimplicit3(@(xq,yq,zq) (interp3((1:S(1))/S(1),(1:S(2))/S(2),(1:S(3))/S(3), bigarray-level*maxamplitude, xq, yq, zq)),'EdgeColor','none','FaceAlpha',.2);
            hold on;
        end
        hold off;
        %    isosurface((1:S(1))/S(1),(1:S(2))/S(2),(1:S(3))/S(3), bigarray, 0.01);
        pause();
    end
end
end