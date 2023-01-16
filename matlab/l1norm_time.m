left_dim = 128;
right_dim = 256;

Al = read_bin('../data/test3curl_left.dat', left_dim / 2 + 1, left_dim + 1);
Ar = read_bin('../data/test3curl_right.dat', right_dim / 2 + 1, right_dim + 1);
Al1 = zeros(1,length(Al(1,1,:)));

Al_left = read_bin('../data/test3lcurl_left.dat', left_dim / 2 + 1, left_dim + 1);
Ar_left = read_bin('../data/test3lcurl_right.dat', left_dim / 2 + 1, left_dim + 1);
Al1_left = zeros(1,length(Al_left(1,1,:)));

Al_right = read_bin('../data/test3rcurl_left.dat', right_dim / 2 + 1, right_dim + 1);
Ar_right = read_bin('../data/test3rcurl_right.dat', right_dim / 2 + 1, right_dim + 1);
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
legend('left scale', 'right scale', 'multiscale')