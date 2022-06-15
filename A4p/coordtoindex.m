function ms = coordtoindex(mi, S)
ms = mi(:,1)+mi(:,2)*S(1)+mi(:,3)*S(1)*S(2);
end