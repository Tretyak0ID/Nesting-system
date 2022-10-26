A = read_bin('../data/test1.b', 129, 129);

axis([0 129 0 129 0 10000])
for i = 1:length(A(1,1,:))
    zlim([9600 11000])
    hold on; grid on
    surf(X, Y, A(:, :, i))
    pause(0.1)
    if (i < length(A(1,1,:)))
        cla;
    end
end