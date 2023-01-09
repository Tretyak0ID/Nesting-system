left_dim = 128;
right_dim = 64;
Al = read_bin('../data/test5_128x64_dynamic_curl_left.dat', left_dim / 2 + 1, left_dim + 1);
Ar = read_bin('../data/test5_128x64_dynamic_curl_right.dat', right_dim / 2 + 1, right_dim + 1);

x_fix_l = floor(length(Al(:, 1, 1)) / 2) + 1;
x_fix_r = floor(length(Ar(:, 1, 1)) / 2) + 1;
y_fix = floor(length(Al(1, :, 1)) / 2);
time_step = 25;

xl = [0 : right_dim / left_dim : right_dim / 2];
xr = [right_dim / 2 : 1 : right_dim];
yl = [0 : right_dim / left_dim : right_dim];
yr = [0 : 1 : right_dim];

%field
hold on; grid on;

% plot(xl, Al(:, y_fix, time_step))
% plot(xr, Ar(:, y_fix, time_step))

plot(yl, Al(x_fix_l, :, time_step) - Al(x_fix, :, 1))
plot(yr, Ar(x_fix_r, :, time_step) - Ar(x_fix, :, 1))

legend('128', '64')
%{
figure;
for i = 1:1:length(Al(1,1,:))
     %zlim([9000 11000])  %comment for 2d field
    plot(yl, Al(x_fix, :, i))
    hold on; grid on;
    plot(yr, Ar(x_fix, :, i))
    pause(0.1)
    if (i < length(Al(1,1,:)))
        cla;
    end
end
%}