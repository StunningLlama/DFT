H = sym([0 0 -4*sqrt(6) 0; 0 0 0 0; -4*sqrt(6) 0 0 0; 0 0 0 0])
P = sym([1/sqrt(2) 0 1/sqrt(2) 0; 0 1 0 0; 1/sqrt(2) 0 -1/sqrt(2) 0; 0 0 0 1])
inv(P)*H*P
P*H*inv(P)