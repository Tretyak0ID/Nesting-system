dim1 = 192;
dim2 = 192 * 2;
dim22 = 192 * 2;

A11 = read_bin('../data/test54_11.dat', dim1 / 3 + 1, dim1 / 2 + 1);
A12 = read_bin('../data/test54_12.dat', dim22 / 3 + 1, dim22 / 2 + 1);
A21 = read_bin('../data/test54_21.dat', dim1 / 3 + 1, dim1 / 2 + 1);
A22 = read_bin('../data/test54_22.dat', dim1 / 3 + 1, dim1 / 2 + 1);
A31 = read_bin('../data/test54_31.dat', dim1 / 3 + 1, dim1 / 2 + 1);
A32 = read_bin('../data/test54_32.dat', dim1 / 3 + 1, dim1 / 2 + 1);

axis([0 dim1 0 dim1 0 10000])
x1 = [0 : 1 : dim1 / 3];
x2 = [dim1 / 3 : 1 : 2 * dim1 / 3];
x3 = [2 * dim1 / 3 : 1 : dim1];
y1 = [0 : 1 : dim1 / 2];
y2 = [dim1 / 2 : 1 : dim1];

x11 = [0 : dim1 / dim22 : dim1 / 3];
y11 = [0 : dim1 / dim22 : dim1 / 2];
x22 = [dim1 / 3 : dim1 / dim22 : 2 * dim1 / 3];
y22 = [dim1 / 2 : dim1 / dim22 : dim1];
x33 = [2 * dim1 / 3 : dim1 / dim22 : dim1];

[X11, Y11] = meshgrid(x1, y1);
[X12, Y12] = meshgrid(x11, y22);
[X21, Y21] = meshgrid(x2, y1);
[X22, Y22] = meshgrid(x2, y2);
[X31, Y31] = meshgrid(x3, y1);
[X32, Y32] = meshgrid(x3, y2);

for i = 1:1:length(A11(1,1,:))
    %zlim([min(min(min(A32(:, :, :)))) max(max(max(A31(:,:,:))))])  %for 3d field
    %zlim([9000 11000])  %for 3d field
    surf(X11, Y11, A11(:, :, i)')
    title(i)
    hold on; grid on;
    surf(X12, Y12, A12(:, :, i)')
    surf(X21, Y21, A21(:, :, i)')
    surf(X22, Y22, A22(:, :, i)')
    surf(X31, Y31, A31(:, :, i)')
    surf(X32, Y32, A32(:, :, i)')
    colormap turbo;
    pause(0.1)
    if (i < length(A11(1,1,:)))
        cla;
    end
end