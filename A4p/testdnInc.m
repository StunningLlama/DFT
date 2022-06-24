BB = [];
global gbl_activeindex;
global gbl_activelength;
kS = [4; 4; 1];

for row = [1:prod(kS)]
    BB = [BB; norm(B([gbl_activeindex(row): gbl_activeindex(row)+gbl_activelength(row)-1],1))];
end


kms=[0:prod(kS)-1]';
km1=rem(kms,kS(1));
km2=rem(floor(kms/kS(1)),kS(2));
km3=rem(floor(kms/(kS(1)*kS(2))),kS(3));
kM=[km1, km2, km3];