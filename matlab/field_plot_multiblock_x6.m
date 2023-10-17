close all;
clear all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
%--------------------------------------------------------------------------


dim1 = 128;
dim22 = dim1 * 2;
L_max = 40000;
test_name = 'test5_x4';

A11 = read_bin('../data/test51_11.dat', dim1 / 2 + 1, dim1 / 2 + 1);
A12 = read_bin('../data/test51_12.dat', dim22 / 2 + 1, dim22 / 2 + 1);
A21 = read_bin('../data/test51_21.dat', dim1 / 2 + 1, dim1 / 2 + 1);
A22 = read_bin('../data/test51_22.dat', dim1 / 2 + 1, dim1 / 2 + 1);
%A31 = read_bin('../data/test5_31.dat', dim1 / 3 + 1, dim1 / 2 + 1);
%A32 = read_bin('../data/test5_32.dat', dim1 / 3 + 1, dim1 / 2 + 1);

axis([0 dim1 0 dim1 0 10000])
x1 = [0 : L_max / dim1 : L_max / 2];
x2 = [L_max / 2 : L_max/dim1 : L_max];
y1 = [0 : L_max / dim1 : L_max / 2];
y2 = [L_max / 2 : L_max / dim1 : L_max];

x11 = [0 : L_max / dim1 / 2 : L_max / 2];
y11 = [0 : L_max / dim1 / 2 : L_max / 2];
x22 = [L_max / 2 : L_max / dim1 / 2 : L_max];
y22 = [L_max / 2 : L_max / dim1 / 2 : L_max];

[X11, Y11] = meshgrid(x1, y1);
[X12, Y12] = meshgrid(x11, y22);
[X21, Y21] = meshgrid(x2, y1);
[X22, Y22] = meshgrid(x2, y2);
%[X31, Y31] = meshgrid(x3, y1);
%[X32, Y32] = meshgrid(x3, y2);

% for i = 1:1:length(A11(1,1,:))
%     %zlim([min(min(min(A32(:, :, :)))) max(max(max(A31(:,:,:))))])  %for 3d field
%     %zlim([9000 11000])  %for 3d field
%     surf(X11, Y11, A11(:, :, i)')
%     hold on; grid on;
%     surf(X12, Y12, A12(:, :, i)')
%     surf(X21, Y21, A21(:, :, i)')
%     surf(X22, Y22, A22(:, :, i)')
%     colormap turbo;
%     pause(0.1)
%     if (i < length(A11(1,1,:)))
%         cla;
%     end
% end
hold on; grid on;
caxis([-5e-5 5e-5]);
c = colorbar;
xlim([0 L_max])
ylim([0 L_max])
pcolor(X11, Y11, A11(:, :, 173)')
%caxis([-5e-5 5e-5]);
pcolor(X12, Y12, A12(:, :, 173)')
%caxis([-5e-5 5e-5]);
pcolor(X21, Y21, A21(:, :, 173)')
%caxis([-5e-5 5e-5]);
pcolor(X22, Y22, A22(:, :, 173)')
%caxis([-5e-5 5e-5]);
shading flat;
xlabel('x, km')
ylabel('y, km')
title("t = 15 days, $\Delta x_{max}$ = 312 km")

% cla;
% hold on; grid on;
% c = colorbar;
% c.Limits = [-5e-5 5e-5];
% 
% contourf(X11, Y11, A11(:, :, 144)', 8)
% caxis([-5e-5 5e-5]);
% contourf(X12, Y12, A12(:, :, 144)', 8)
% caxis([-5e-5 5e-5]);
% contourf(X21, Y21, A21(:, :, 144)', 8)
% caxis([-5e-5 5e-5]);
% contourf(X22, Y22, A22(:, :, 144)', 8)
% caxis([-5e-5 5e-5]);
% colormap jet;
% shading flat;
% save_file_name = strcat('../data/', test_name, '_12.png');
% saveas(gcf, save_file_name);
% 
% cla;
% hold on; grid on;
% c = colorbar;
% c.Limits = [-5e-5 5e-5];

% contourf(X11, Y11, A11(:, :, 173)', 8)
% caxis([-5e-5 5e-5]);
% contourf(X12, Y12, A12(:, :, 173)', 8)
% caxis([-5e-5 5e-5]);
% contourf(X21, Y21, A21(:, :, 173)', 8)
% caxis([-5e-5 5e-5]);
% contourf(X22, Y22, A22(:, :, 173)', 8)
% caxis([-5e-5 5e-5]);
% colormap jet;
% shading flat;
% save_file_name = strcat('../data/', test_name, '_15.png');
% saveas(gcf, save_file_name);


%--------------------------------------------------------------------------
set(gca, 'FontSize', 20)

exportgraphics(gcf,'Kelhelmref12.pdf')