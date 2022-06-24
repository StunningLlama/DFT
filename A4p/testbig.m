X=[8 8 8; 8+2 8 8];
Z = [1 1];

setup(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), [2; 2; 2], false);
Wsmall=initializeRandomState();
%Wsmall = iterate(20);
getE(Wsmall)

setupLargecell(X, 1, [1 1], [48; 48; 48], diag([16 16 16]), [2; 2; 2], false);
Wbig = initializeZeroState();
for k = [1:8]
    Wbig{1} = setLargeCell(Wbig{1}, k, k, Wsmall{k});
end
%W = iterate(40);
getE(Wbig)





% R = diag([16 16 16]);
% kS = [2; 2; 2];
% kms=[0:prod(kS)-1]';
% km1=rem(kms,kS(1));
% km2=rem(floor(kms/kS(1)),kS(2));
% km3=rem(floor(kms/(kS(1)*kS(2))),kS(3));
% kM=[km1, km2, km3];
% 
% Xsmall = X;
% Zsmall = Z;
% X = [];
% Z = [];
% for k = [1:prod(kS)]
%     offsets = kM*(R');
%     X = [X; Xsmall + offsets(k,:)];
%     Z = [Z Zsmall];
% end
% X
% Z
% 
% setup(X, 1*8, Z, [48; 48; 48]*2, diag([16 16 16])*2, [1; 1; 1], false);
% W = iterate(40);
% getE(W)