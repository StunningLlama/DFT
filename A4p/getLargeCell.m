function data = getLargeCell(W, row, col)
global gbl_Ns;
global gbl_kS;
Ns_small = gbl_Ns/gbl_kS;
colrange = (col-1)*Ns_small + [1:Ns_small];
rowrange = [activeindex(row): activeindex(row)+activelength(row)];
data = W(rowrange,colrange);
end