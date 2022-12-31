dim = 128;
A = read_bin('../data/test5_128x128_curl_dynamic.dat', dim + 1, dim + 1);

x_fix = floor(length(Al(:, 1, 1)) / 2);
y_fix = floor(length(Al(1, :, 1)) / 2);
time_step = 2;

xl = [0 : 1 : dim];
yl = [0 : 1 : dim];

%field
hold on; grid on;

%plot(xl, Al(:, y_fix, time_step))
%plot(xr, Ar(:, y_fix, time_step))

plot(y, A(x_fix, :, time_step) - A(x_fix, :, 1))