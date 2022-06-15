B = randn(4);
B = B*B';
dB = randn(4);
dB = 0.0001*(dB*dB');
in = randn(4)+i*randn(4);
Q(in, B+dB) - Q(in, B)
dQ(in, B, dB)