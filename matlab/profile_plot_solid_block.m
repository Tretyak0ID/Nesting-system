dim = 128;
A = read_bin('../data/test2h.dat', dim + 1, dim + 1);

x_fix = floor(length(A(:, 1, 1)) / 2);
y_fix = floor(length(A(1, :, 1)) / 2);
time_step = 2;

x = [0 : 1 : dim];
y = [0 : 1 : dim];

%field
hold on; grid on;

%plot(xl, Al(:, y_fix, time_step))
%plot(xr, Ar(:, y_fix, time_step))

plot(y, A(x_fix, :, end))