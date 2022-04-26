syms P Q P0 Q0 t C
C = Q0-(log(P0^2+1))/2
P = exp(-2*C)*t+P0
Q = log((exp(-2*C)*t+P0)^2+1)/2+C

A = diff(Q,t)-P*exp(-2*Q)
B = diff(P,t)-(P^2+1)*exp(-2*Q)
