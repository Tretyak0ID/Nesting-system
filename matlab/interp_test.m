I1 = read_bin('../data/test_interp1.dat', 128 + 1, 0 + 1);
I2 = read_bin('../data/test_interp2.dat', 64 + 1, 0 + 1);

hold on; grid on;
x1 = [0:0.5:64];
x2 = [0:1:64];
plot(x1, I1)
plot(x2, I2)

I3 = read_bin('../data/test_interp3.dat', 64 + 1, 0 + 1);
I4 = read_bin('../data/test_interp4.dat', 128 + 1, 0 + 1);

figure
hold on; grid on;
plot(x2, I3)
plot(x1, I4)