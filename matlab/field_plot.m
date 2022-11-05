A = read_bin('../data/test_interp.dat', 65, 1);

plot(A)

figure;
axis([0 129 0 129 0 10000])
for i = 1:5:length(A(1,1,:))
    %zlim([9600 11000])  %comment for 2d field
    hold on; grid on
    surf(A(:, :, i))
    pause(0.1)
    if (i < length(A(1,1,:)))
        cla;
    end
end