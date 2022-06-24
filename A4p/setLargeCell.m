function out = setLargeCell(W, row, col, data)
global gbl_Ns;
global gbl_kS;
global gbl_activeindex;
global gbl_activelength;
global gbl_Ns_small;
out = W;
Ns_small = gbl_Ns_small;
colrange = (col-1)*Ns_small + [1:Ns_small];
rowrange = [gbl_activeindex(row): gbl_activeindex(row)+gbl_activelength(row)-1];
out(rowrange,colrange) = data;
end