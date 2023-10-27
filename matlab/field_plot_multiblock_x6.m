%--------------------------------------------------------------------------


dim1 = 128;
dim22 = dim1 * 2;
L_max = 40000;

A11 = read_bin('../data/test51_11.dat', dim1 / 2 + 1, dim1 / 2 + 1);
A12 = read_bin('../data/test51_12.dat', dim22 / 2 + 1, dim22 / 2 + 1);
A21 = read_bin('../data/test51_21.dat', dim1 / 2 + 1, dim1 / 2 + 1);
A22 = read_bin('../data/test51_22.dat', dim1 / 2 + 1, dim1 / 2 + 1);
%A31 = read_bin('../data/test5_31.dat', dim1 / 3 + 1, dim1 / 2 + 1);
%A32 = read_bin('../data/test5_32.dat', dim1 / 3 + 1, dim1 / 2 + 1);

axis([0 L_max 0 L_max 0 10000])
x1 = [0 : L_max / dim1 : L_max / 2];
x2 = [L_max / 2 : L_max/dim1 : L_max];
y1 = [0 : L_max / dim1 : L_max / 2];
y2 = [L_max / 2 : L_max / dim1 : L_max];

x11 = [0 : L_max / dim1 / (dim22 / dim1) : L_max / 2];
y11 = [0 : L_max / dim1 / (dim22 / dim1) : L_max / 2];
x22 = [L_max / 2 : L_max / dim1 / (dim22 / dim1) : L_max];
y22 = [L_max / 2 : L_max / dim1 / (dim22 / dim1) : L_max];

[X11, Y11] = meshgrid(x1, y1);
[X12, Y12] = meshgrid(x11, y22);
[X21, Y21] = meshgrid(x2, y1);
[X22, Y22] = meshgrid(x2, y2);
%[X31, Y31] = meshgrid(x3, y1);
%[X32, Y32] = meshgrid(x3, y2);

for i = 1:1:length(A11(1,1,:))
    %zlim([min(min(min(A32(:, :, :)))) max(max(max(A31(:,:,:))))])  %for 3d field
    %zlim([9000 11000])  %for 3d field
    hold on; grid on;
    c = colorbar;
    caxis([-5e-5 5e-5]);
    pcolor(X11, Y11, A11(:, :, i)')
    pcolor(X12, Y12, A12(:, :, i)')
    pcolor(X21, Y21, A21(:, :, i)')
    pcolor(X22, Y22, A22(:, :, i)')
    shading interp
    pause(0.1)
    if (i < length(A11(1,1,:)))
        cla;
    end
end
% hold on; grid on;
% caxis([-5e-5 5e-5]);
% c = colorbar;
% xlim([0 L_max])
% ylim([0 L_max])
% pcolor(X11, Y11, A11(:, :, 5)')
% pcolor(X12, Y12, A12(:, :, 5)')
% pcolor(X21, Y21, A21(:, :, 5)')
% pcolor(X22, Y22, A22(:, :, 5)')
% shading flat;
% xlabel('x, km')
% ylabel('y, km')
% title("t = 15 days, $\Delta x_{max}$ = 312 km")


