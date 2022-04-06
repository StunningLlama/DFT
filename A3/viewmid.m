% Function to view slices of three dimensional data sets
%
% Usage: view(dat,S)
%
% dat: 3d data set (any shape) of total size prod(S)=S(1)*S(2)*S(3)
% S: dimensions of dat in a 3-vector
function viewmid(dat,S)
fprintf('\nRemember to hit <enter> or <spacebar> after each plot!\n\n');
for k=1:3
text=sprintf("m%d=%d slice",k,floor(S(k)/2));
fprintf("%s...\n",text);
title(text);
if k==1
xlabel("m3 ->"); ylabel("m2 ->");
elseif k==2
xlabel("m3 ->"); ylabel("m1 ->");
elseif k==3
xlabel("m2 ->"); ylabel("m1 ->");
end
mesh(slice(dat,S,S(k)/2,k)); pause;
end