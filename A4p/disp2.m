function disp2(A)
vertdim = size(A,1);
horizdim = size(A,2);
disp(strcat(inputname(1), ": ", num2str(vertdim), "*", num2str(horizdim), " matrix."));
maxdim = 6;
if (horizdim < maxdim)
    horizlist = [1:horizdim];
else
    horizlist = [1:maxdim/2 horizdim-(maxdim/2-1):horizdim];
end
if (vertdim < maxdim)
    vertlist = [1:vertdim];
else
    vertlist = [1:maxdim/2 vertdim-(maxdim/2-1):vertdim];
end
B = A(vertlist,horizlist);
disp(B);
end