% Function to view slices of three dimensional data sets
%
% Usage: view(dat,S)
%
% dat: 3d data set (any shape) of total size prod(S)=S(1)*S(2)*S(3)
% S: dimensions of dat in a 3-vector
function viewdft(dat,S)
fprintf('\nRemember to hit <enter> or <spacebar> after each plot!\n\n');
for k=1:3
if k==1
fprintf('m1=0 slice...\n');
title("m1=0 slice"); xlabel("m3 ->"); ylabel("m2 ->");
elseif k==2
fprintf('m2=0 slice...\n');
title("m2=0 slice"); xlabel("m3 ->"); ylabel("m1 ->");
elseif k==3
fprintf('m3=0 slice...\n');
title("m3=0 slice"); xlabel("m2 ->"); ylabel("m1 ->");
end
mesh(slice(dat,S,1,k)); pause;
end
end