Al = read_bin('../data/test1h_left.dat', 129, 257);
Ar = read_bin('../data/test1h_right.dat', 65, 129);

axis([0 128 0 128 0 10000])
xl = [0 : 0.5 : 64];
xr = [64 : 1 : 128];
yl = [0 : 0.5 : 128];
yr = [0 : 1 : 128];

[Xl, Yl] = meshgrid(yl, xl);
[Xr, Yr] = meshgrid(yr, xr);

for i = 1:1:length(Al(1,1,:))
    zlim([9600 11000])  %comment for 2d field
    hold on; grid on
    surf(Xl, Yl, Al(:, :, i))
    surf(Xr, Yr, Ar(:, :, i))
    pause(0.1)
    if (i < (length(Al(1,1,:)) - 6))
        cla;
    end
end