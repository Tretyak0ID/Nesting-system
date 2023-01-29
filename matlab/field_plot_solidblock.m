dim = 192;
A = read_bin('../data/test6_192_curl.dat', dim + 1, dim + 1);

axis([0 dim 0 dim 0 10000])
x = [0 : 1 : dim];
y = [0 : 1 : dim];

[X, Y] = meshgrid(x, y);

%field
for i = 1:1:length(A(1, 1, :))
    zlim([min(min(min(A))) max(max(max(A)))])  %comment for 2d field
    surf(X, Y, A(:, :, i)')
    hold on; grid on;
    shading flat
    colormap turbo
    pause(0.2)
    if (i < length(A(1, 1, :)))
        cla;
    end
end

%difference
%{
figure;
for i = 1:1:length(A(1, 1, :))
    %zlim([9000 11000])  %comment for 2d field
    surf(X, Y, A(:, :, i) -  A(:, :, 1))
    hold on; grid on;
    shading flat
    pause(0.1)
    if (i < length(A(1, 1, :)))
        cla;
    end
end
%}