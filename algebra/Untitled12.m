syms Z W E
E = Z*((Z^2*(16*Z-22)+W*Z*(50*Z-88)+W^2*(64*Z-109))/(16*Z^2+70*W*Z+96*W^2));
Ew = diff(E,W);
Ez = diff(E,Z);
sol = solve(Ew==0, Ez==0);
for i = [1:size(sol.W, 1)]
Emin = subs(E, [W Z], [sol.W(i) sol.Z(i)]);
[double(sol.W(i)) double(sol.Z(i)) double(Emin)]
end