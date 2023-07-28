dim = 192 * 2;
L_max = 40000;
A = read_bin('../data/test3curl.dat', dim + 1, dim + 1);
test_name = 'test5386x386';

axis([0 dim 0 dim 0 10000])
x = [0 : L_max / dim : L_max];
y = [0 : L_max / dim : L_max];

[X, Y] = meshgrid(x, y);

%field
% hold on; grid on;
% colorbar;
%  for i = 1:1:length(A(1, 1, :))
%      zlim([min(min(min(A))) max(max(max(A)))])  %comment for 2d field
%      surf(X, Y, A(:, :, i)')
%      shading flat
%      colormap turbo
%      pause(0.2)
%      if (i < length(A(1, 1, :)))
%          cla;
%      end
%  end

hold on; grid on;
%caxis([-2e-5 4.5e-5]);
%c = colorbar;
%zlim([min(min(min(A22(:, :, :)))) max(max(max(A22(:,:,:))))])  %for 3d field
%c.Limits = [-2e-5 4.5e-5];
contourf(X, Y, A(:, :, 1)', 15)
colormap jet
xlabel('x, км')
ylabel('y, км')