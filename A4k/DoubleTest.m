setupSmallGe(24);
ewald()
W = iterate(20);
E0 = getE(W)

global gbl_X;
%visualize(W, gbl_X);
global gbl_S;
global gbl_M;
global gbl_kpoints;
global gbl_kvectors;
global gbl_kR;
global gbl_Ns;
IWnew = zeros(prod(gbl_S*2), gbl_kpoints*gbl_Ns);
indices = coordtoindex(gbl_M, gbl_S);
coords = indextocoord(indices, gbl_S);
for k = [1:gbl_kpoints]
    for j = [1:gbl_Ns]
        kvec = gbl_kR(k,:);
        IW = cI(W{k}(:,j), k);
        (k-1)*gbl_Ns+(j-1)+1
        for dx = [0:1]
            for dy = [0:1]
                for dz = [0:1]
                    newcoords = coords + gbl_S'.*[dx dy dz];
                    %disp2(newcoords);
                    newindices = coordtoindex(newcoords, gbl_S*2)+1;
                    %disp2(newindices);
                    karray = ones(size(newcoords,1),1)*kvec;
                    coordarray = newcoords./(ones(size(newcoords,1),1)*gbl_S');
                    % disp2(karray);
                    % disp2(coordarray);
%                     pause();
                    phases = exp(2*pi*i*sum(coordarray.*karray,2));
                    %disp2(phases);
                    %pause();
                    IWnew(newindices,(k-1)*gbl_Ns+(j-1)+1) = IW.*phases;
                end
            end
        end
    end
end
setupBigGe(24);
Wnew = {cJcomp(IWnew, 1)};
global gbl_X;
%visualize(Wnew, gbl_X);
disp("Energy");
E1 = getE(Wnew)
ratio = E1/E0
%Wnew=sd(Wnew,5);
%Wnew = orthonormalize(Wnew); %# Restart as orthonormal functions
Wnew=pccg(Wnew,50,1); %# 50 iterations of pclm from same W
E2 = getE(Wnew)
ratio2 = E2/E0
