close all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
%--------------------------------------------------------------------------


dim1 = 96;
dim2 = 96 * 2;
dim22 = 96 * 2;
L_max = 4000;
%dim3 = 192;
%A1 = read_bin('../data/test1h_1.dat', dim1 + 1, dim1 / 2 + 1);
%A2 = read_bin('../data/test1_2.dat', dim2 + 1, dim2 / 2 + 1);
A11 = read_bin('../data/test6_96_11_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A12 = read_bin('../data/test6_96_12_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A13 = read_bin('../data/test6_96_13_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A21 = read_bin('../data/test6_96_21_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A22 = read_bin('../data/test6_96_22_curl.dat', dim22 / 3 + 1, dim22 / 3 + 1);
A23 = read_bin('../data/test6_96_23_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A31 = read_bin('../data/test6_96_31_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A32 = read_bin('../data/test6_96_32_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A33 = read_bin('../data/test6_96_33_curl.dat', dim1 / 3 + 1, dim1 / 3 + 1);

axis([0 L_max 0 L_max 0 10000])
x1 = [0 : L_max / dim1 : L_max / 3];
x2 = [L_max / 3 : L_max / dim1 : 2 * L_max / 3];
x3 = [2 * L_max / 3 : L_max / dim1 : L_max];
y1 = [0 : L_max / dim1 : L_max / 3];
y2 = [L_max / 3 : L_max / dim1 : 2 * L_max / 3];
y3 = [2 * L_max / 3 : L_max/ dim1 : L_max];

x11 = [0 : L_max / dim1 / 2 : L_max / 3];
y11 = [0 : L_max / dim1 / 2 : L_max / 3];
x22 = [L_max / 3 : L_max / dim1 / 2 : 2 * L_max / 3];
y22 = [L_max / 3 : L_max / dim1 / 2 : 2 * L_max / 3];
x33 = [2 * L_max / 3 : L_max / dim1 / 2 : L_max];
y33 = [2 * L_max / 3 : L_max / dim1 / 2 : L_max];

[X11, Y11] = meshgrid(x1, y1);
[X12, Y12] = meshgrid(x1, y2);
[X13, Y13] = meshgrid(x1, y3);
[X21, Y21] = meshgrid(x2, y1);
[X22, Y22] = meshgrid(x22, y22);
[X23, Y23] = meshgrid(x2, y3);
[X31, Y31] = meshgrid(x3, y1);
[X32, Y32] = meshgrid(x3, y2);
[X33, Y33] = meshgrid(x3, y3);

xx1 = [0 : 1 : dim1];
xx2 = [0 : 1/2 : dim1];
yy1 = [0 : 1 : dim1 / 2];
yy2 = [dim1 / 2 : 1/2 : dim1];
[X1, Y1] = meshgrid(xx1, yy1);
[X2, Y2] = meshgrid(xx2, yy2);

% for i = 1:1:length(A11(1,1,:))
%     zlim([min(min(min(A22(:, :, :)))) max(max(max(A22(:,:,:))))])  %for 3d field
%     %zlim([9000 11000])  %for 3d field
%     %surf(X1, Y1, A1(:, :, i)')
% 
%     surf(X11, Y11, A11(:, :, i)')
%     title(i)
%     hold on; grid on;
%     colorbar;
%     %surf(X2, Y2, A2(:, :, i)')
%     surf(X12, Y12, A12(:, :, i)')
%     surf(X13, Y13, A13(:, :, i)')
%     surf(X21, Y21, A21(:, :, i)')
%     surf(X22, Y22, A22(:, :, i)')
%     surf(X23, Y23, A23(:, :, i)')
%     surf(X31, Y31, A31(:, :, i)')
%     surf(X32, Y32, A32(:, :, i)')
%     surf(X33, Y33, A33(:, :, i)')
%     colormap jet;
%     pause(0.1)
%     if (i < length(A11(1,1,:)))
%         cla;
%     end
% end

hold on; grid on;
c = colorbar;
%zlim([min(min(min(A22(:, :, :)))) max(max(max(A22(:,:,:))))])  %for 3d field
%c.Limits = [9200 10100];
pcolor(X11, Y11, A11(:, :, 75)')
pcolor(X12, Y12, A12(:, :, 75)')
pcolor(X13, Y13, A13(:, :, 75)')
pcolor(X21, Y21, A21(:, :, 75)')
pcolor(X22, Y22, A22(:, :, 75)')
pcolor(X23, Y23, A23(:, :, 75)')
pcolor(X31, Y31, A31(:, :, 75)')
pcolor(X32, Y32, A32(:, :, 75)')
pcolor(X33, Y33, A33(:, :, 75)')
%caxis([9200 10100]);
shading interp;
xlabel('x, km')
ylabel('y, km')
title("t = 10 hours, $\Delta x_{max}$ = 42 km")


%--------------------------------------------------------------------------
set(gca, 'FontSize', 20)