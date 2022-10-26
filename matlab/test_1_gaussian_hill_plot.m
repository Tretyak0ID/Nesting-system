A = read_bin('../data/test1.b', 129, 129);

hold on; grid on;
for i = 1:length(A(1,1,:))
    pause(0.001)
    surf(A(:, :, i))
    drawnow
end