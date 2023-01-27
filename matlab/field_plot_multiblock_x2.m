dim1 = 96;
dim2 = 96;

A1 = read_bin('../data/test6_1_curl.dat', dim1 / 2 + 1, dim1 + 1);
A2 = read_bin('../data/test6_2_curl.dat', dim2 / 2 + 1, dim2 + 1);

axis([0 dim1 0 dim1 0 10000])
x1 = [0 : 1 : dim1 / 2];
x2 = [dim1 / 2 : dim1 / dim2 : dim1];
y1 = [0 : 1 : dim1];
y2 = [0 : dim1 / dim2 : dim1];

[X1, Y1] = meshgrid(x1, y1);
[X2, Y2] = meshgrid(x2, y2);

for i = 1:1:length(A1(1,1,:))
    %zlim([min(min(min(A32(:, :, :)))) max(max(max(A31(:,:,:))))])  %for 3d field
    %zlim([9000 11000])  %for 3d field
    surf(X1, Y1, A1(:, :, i)')
    title(i)
    hold on; grid on;
    surf(X2, Y2, A2(:, :, i)')
    colormap turbo;
    shading interp;
    pause(0.1)
    if (i < length(A1(1,1,:)))
        cla;
    end
end