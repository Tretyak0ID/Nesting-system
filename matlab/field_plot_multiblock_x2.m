dim1 = 256;
dim2 = 128;
L_max = 40000;

A1 = read_bin('../data/diffusion_left.dat', dim1 / 2 + 1, dim1 + 1);
A2 = read_bin('../data/diffusion_right.dat', dim2 / 2 + 1, dim2 + 1);

axis([0 dim1 0 dim1 0 10000])
x1 = [0 : L_max / dim1 : L_max / 2];
x2 = [L_max / 2:  L_max / dim2 : L_max];
y1 = [0 : L_max / dim1 : L_max];
y2 = [0 : L_max / dim2 : L_max];

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
    pause(0.001)
    if (2 * i < length(A1(1,1,:)))
        cla;
    end
end

% surf(X1, Y1, A1(:, :, 1)')
% hold on; grid on;
% surf(X2, Y2, A2(:, :, 1)')
% colormap jet;
% xlabel('x, км')
% ylabel('y, км')
% c = colorbar;
% %c.Limits = [1e4 1.1e4];