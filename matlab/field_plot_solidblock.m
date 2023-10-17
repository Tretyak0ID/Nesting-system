close all;
clear all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
%--------------------------------------------------------------------------

dim = 96;
L_max = 4000;
A = read_bin('../data/test6_96_h.dat', dim + 1, dim + 1);

%axis([0 dim 0 dim 0 10000])
x = [0 : L_max / dim : L_max];
y = [0 : L_max / dim : L_max];

[X, Y] = meshgrid(x, y);

% hold on; grid on;
%  for i = 1:1:length(A(1, 1, :))
%      %zlim([min(min(min(A))) max(max(max(A)))])  %comment for 2d field
%      surf(X, Y, A(:, :, i)' - A(:,:,1)')
%      title(i)
%      shading flat
%      colormap turbo
%      pause(0.2)
%      if (i < length(A(1, 1, :)))
%          cla;
%      end
%  end

hold on; grid on;
caxis([9200 10100]);
c = colorbar;
%zlim([min(min(min(A22(:, :, :)))) max(max(max(A(:,:,:))))])  %for 3d field
c.Limits = [9200 10100];
pcolor(X, Y, A(:, :, 216)') %173 138
shading flat;
colormap jet
xlabel('x, km')
ylabel('y, km')
title("t = 216 hours, $\Delta x$ = 42 km")



%--------------------------------------------------------------------------
set(gca, 'FontSize', 20)

exportgraphics(gcf,'Kelhelmref12.pdf')