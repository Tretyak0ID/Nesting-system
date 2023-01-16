function out = plot_l1norm(multiscale_test, leftsclae_test, rightscale_test, left_dim, right_dim)

    left_multiscale_name = strcat('../data/', multiscale_test, '_left.dat');
    right_multiscale_name = strcat('../data/', multiscale_test, '_right.dat');
    Al = read_bin(left_multiscale_name, left_dim / 2 + 1, left_dim + 1);
    Ar = read_bin(right_multiscale_name, right_dim / 2 + 1, right_dim + 1);
    Al1 = zeros(1,length(Al(1,1,:)));
    
    left_leftscale_name = strcat('../data/', leftsclae_test, '_left.dat');
    right_leftscale_name = strcat('../data/', leftsclae_test, '_right.dat');
    Al_left = read_bin(left_leftscale_name, left_dim / 2 + 1, left_dim + 1);
    Ar_left = read_bin(right_leftscale_name, left_dim / 2 + 1, left_dim + 1);
    Al1_left = zeros(1,length(Al_left(1,1,:)));
    
    left_rightscale_name = strcat('../data/', rightscale_test, '_left.dat');
    right_rightscale_name = strcat('../data/', rightscale_test, '_right.dat');
    Al_right = read_bin(left_rightscale_name, right_dim / 2 + 1, right_dim + 1);
    Ar_right = read_bin(right_rightscale_name, right_dim / 2 + 1, right_dim + 1);
    Al1_right = zeros(1,length(Al_right(1,1,:)));
    
    time = [1 : length(Al1)];
    time_left = [1 : length(Al1_left)];
    time_right = [1 :left_dim / right_dim: length(Al1)];
    
    for i = 1: length(Al1)
        M1 = max(max(abs(Al(:,:,i))));
        M2 = max(max(abs(Ar(:,:,i))));
        Al1(i) = max(M1,M2);
    end
    
    for i = 1: length(Al1_left)
        M1 = max(max(abs(Al_left(:,:,i))));
        M2 = max(max(abs(Ar_left(:,:,i))));
        Al1_left(i) = max(M1,M2);
    end
    
    for i = 1: length(Al1_right)
        M1 = max(max(abs(Al_right(:,:,i))));
        M2 = max(max(abs(Ar_right(:,:,i))));
        Al1_right(i) = max(M1,M2);
    end
    
    figure;
    hold on; grid on;
    plot(time_left(1:end), Al1_left(1:end));
    plot(time_right(1:end), Al1_right(1:end));
    plot(time(1:end), Al1(1:end));
    xlabel('t'); ylabel('l1');
    legend('left scale', 'right scale', 'multiscale')
    save_file_name = strcat('../data/', multiscale_test, '_l1_norm.png');
    saveas(gcf, save_file_name);
end