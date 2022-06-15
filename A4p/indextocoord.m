function m = indextocoord(ms, S)
m = [rem(ms,S(1)) rem(floor(ms/S(1)),S(2)) rem(floor(ms/(S(1)*S(2))),S(3))];
end