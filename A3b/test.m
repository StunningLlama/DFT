
uinv = inv(W'*O(W));
n=getn(cI(Y),gbl_f); %# Charge density
%# Basic slice from ‘‘100’’ edge of cell
sl=reshape(n(1:S(1)*S(2)),S(1),S(2));
%# Make and view image
contourf(sl);
%# 110 slice cutting bonds (assumes cube)
sl=reshape(n(find(M(:,2)==M(:,3))),S(1),S(2));
%# Expand by 2, drop data to restore (approximate) aspect ratio
li=find(rem([1:size(sl,1)],3)~=0); sl=sl(li,:);
%# Make and view image
figure(2);
contourf(sl);