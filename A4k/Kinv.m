function out=Kinv(W)
global gbl_Gc;
global gbl_R;
global gbl_kvectors;
global gbl_weights;
global gbl_kpoints;
out = {};

for k = [1:gbl_kpoints]
    kvec = gbl_kvectors(k,:);
    karray = ones(size(W{k}, 1), 1)*kvec;
    
    %out{k} = -(1/gbl_weights(k))*W{k}./(det(gbl_R)*((sum((gbl_Gc+karray).^2, 2)+1)*ones(1,size(W,2))));
    out{k} = (gbl_weights(k))*W{k}.*((sum((gbl_Gc{k}+karray).^2, 2)+1)*ones(1,size(W{k},2)));
end
end