left_dim = 128;
right_dim = 128;
Al = read_bin('../data/test2h_left.dat', left_dim / 2 + 1, left_dim + 1);
Ar = read_bin('../data/test2h_right.dat', right_dim / 2 + 1, right_dim + 1);

x_fix_l = floor(length(Al(:, 1, 1)) / 2) + 1;
x_fix_r = floor(length(Ar(:, 1, 1)) / 2) + 1;
y_fix = floor(length(Al(1, :, 1)) / 2);
time_step = 2;

if (left_dim > right_dim)
    xl = [0 : right_dim / left_dim : right_dim / 2];
    xr = [right_dim / 2 : 1 : right_dim];
    yl = [0 : right_dim / left_dim : right_dim];
    yr = [0 : 1 : right_dim];
else
    xl = [0 : 1 : left_dim / 2];
    xr = [left_dim / 2 : left_dim / right_dim : left_dim];
    yl = [0 : 1 : left_dim];
    yr = [0 : left_dim / right_dim : left_dim];
end

hold on; grid on;
plot(yl, Al(end, :, time_step))
plot(yr, Ar(1, :, time_step))
plot(yl, Al(end, :, 1))
plot(yr, Ar(1, :, 1))
legend('left', 'right')

figure
hold on; grid on;
plot(yl, Al(end, :, time_step) - Al(end, :, 1))
plot(yr, Ar(1, :, time_step) - Ar(1, :, 1))
legend('left', 'right')

% plot(yl, Al(x_fix_l, :, time_step) - Al(x_fix, :, 1))
% plot(yr, Ar(x_fix_r, :, time_step) - Ar(x_fix, :, 1))

%legend('128', '64')

%{
figure;
for i = 1:1:length(Al(1,1,:))
    %ylim([-10^(-4.7) 10^(-4.7)])  %comment for 2d field
    plot(yl, Al(end, :, 2))
    hold on; grid on;
    title(i)
    plot(yr, Ar(1, :, 2))

    pause(0.2)
    if (i < length(Al(1,1,:)))
        cla;
    end
end
%}