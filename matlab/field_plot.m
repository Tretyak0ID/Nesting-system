left_dim = 512;
right_dim = 256;
Al = read_bin('../data/test5curl_left.dat', left_dim / 2 + 1, left_dim + 1);
Ar = read_bin('../data/test5curl_right.dat', right_dim / 2 + 1, right_dim + 1);

axis([0 right_dim 0 right_dim 0 10000])
xl = [0 : right_dim / left_dim : right_dim / 2];
xr = [right_dim / 2 : 1 : right_dim];
yl = [0 : right_dim / left_dim : right_dim];
yr = [0 : 1 : right_dim];

[Xl, Yl] = meshgrid(yl, xl);
[Xr, Yr] = meshgrid(yr, xr);

for i = 1:1:length(Al(1,1,:))
    %zlim([9000 11000])  %comment for 2d field
    surf(Xl, Yl, Al(:, :, i))
    hold on; grid on;
    surf(Xr, Yr, Ar(:, :, i))
    pause(0.1)
    if (i < (length(Al(1,1,:))))
        cla;
    end
end