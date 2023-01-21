dim = 64;
A = read_bin('../data/sbp_test_y.dat', dim + 1, dim + 1);
B = read_bin('../data/cos_field_y.dat', dim + 1, dim + 1);

x_fix = 64;
y_fix = 1;
time_step = 2;

x = [0 : 1 : dim];
y = [0 : 1 : dim];

%field
hold on; grid on;

max(abs(A(:, y_fix) - B(:, y_fix)))