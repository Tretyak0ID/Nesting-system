function out = plotSBPx2multiblok(test_name, left_dim, right_dim, plot3d)
%The function of rendering a field on a two-block SBP-SAT grid

left_block_name = strcat('../data/', test_name, '_left.dat');
right_block_name = strcat('../data/', test_name, '_right.dat');

Al = read_bin(left_block_name, left_dim / 2 + 1, left_dim + 1);
Ar = read_bin(right_block_name, right_dim / 2 + 1, right_dim + 1);

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

figure;
if plot3d == 1
    zlim([min(min(min(Al))) max(max(max(Al)))])  %for 3d field
    surf(Xl, Yl, Al(:, :, 1))
    title(1)
    hold on; grid on;
    surf(Xr, Yr, Ar(:, :, 1))
    colormap turbo
    colorbar;
else
    contourf(Xl, Yl, Al(:, :, 1), 8)
    title(1)
    hold on; grid on;
    contourf(Xr, Yr, Ar(:, :, 1), 8)
    shading flat
    colormap turbo
    colorbar;
end
save_file_name = strcat('../data/', test_name, '_initial_cond.png');
saveas(gcf, save_file_name);
cla;

figure;
if plot3d == 1
    zlim([min(min(min(Al))) max(max(max(Al)))])  %for 3d field
    surf(Xl, Yl, Al(:, :, end))
    title(length(Al(1,1,:)))
    hold on; grid on;
    surf(Xr, Yr, Ar(:, :, end))
    colormap turbo
    colorbar;
else
    contourf(Xl, Yl, Al(:, :, end), 8)
    title(length(Al(1,1,:)))
    hold on; grid on;
    contourf(Xr, Yr, Ar(:, :, end), 8)
    shading flat
    colormap turbo
    colorbar;
end
save_file_name = strcat('../data/', test_name, '_end_time.png');
saveas(gcf, save_file_name);
cla;

figure;
if plot3d == 1
    zlim([min(min(min(Al))) max(max(max(Al)))])  %for 3d field
    surf(Xl, Yl, Al(:, :, floor(length(Al(1,1,:)) / 2)))
    title(floor(length(Al(1,1,:)) / 2))
    hold on; grid on;
    surf(Xr, Yr, Ar(:, :, floor(length(Al(1,1,:)) / 2)))
    colormap turbo
    colorbar;
else
    contourf(Xl, Yl, Al(:, :, floor(length(Al(1,1,:)) / 2)), 8)
    title(floor(length(Al(1,1,:)) / 2))
    hold on; grid on;
    contourf(Xr, Yr, Ar(:, :, floor(length(Al(1,1,:)) / 2)), 8)
    shading flat
    colormap turbo
    colorbar;
end
save_file_name = strcat('../data/', test_name, '_ middle_time.png');
saveas(gcf, save_file_name);


end