left_dim = 128;
right_dim = 256;
Al = read_bin('../data/test2h_left.dat', left_dim / 2 + 1, left_dim + 1);
Ar = read_bin('../data/test2h_right.dat', right_dim / 2 + 1, right_dim + 1);

axis([0 right_dim 0 right_dim 0 10000])
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

[Xl, Yl] = meshgrid(yl, xl);
[Xr, Yr] = meshgrid(yr, xr);

for i = 1:1:length(Al(1,1,:))
    zlim([min(min(min(Al))) max(max(max(Al)))])  %for 3d field
    surf(Xl, Yl, Al(:, :, i))
    title(i)
    hold on; grid on;
    surf(Xr, Yr, Ar(:, :, i))
    shading flat
    pause(0.1)
    if (i < length(Al(1,1,:)))
        cla;
    end
end


%{
for i = 1:1:length(Al(1, 1, :))
    %zlim([0 10 ** (-5)])  %comment for 2d field
    surf(Xl, Yl, Al(:, :, i) -  Al(:, :, 1))
    hold on; grid on;
    surf(Xr, Yr, Ar(:, :, i) - Ar(:, :, 1))
    pause(0.1)
    if (i < length(Al(1, 1, :)))
        cla;
    end
end
%}

%{
surf(Xl, Yl, Al(:, :, 20) - Al(:, :, 1))
hold on; grid on;
surf(Xr, Yr, Ar(:, :, 20) - Ar(:, :, 1))
%}