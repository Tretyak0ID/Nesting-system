dim1 = 192;
dim2 = 192 * 2;

A11 = read_bin('../data/test6_192_11_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A12 = read_bin('../data/test6_192_12_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A13 = read_bin('../data/test6_192_13_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A21 = read_bin('../data/test6_192_21_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A22 = read_bin('../data/test6_192_22_h.dat', dim2 / 3 + 1, dim2 / 3 + 1);
A23 = read_bin('../data/test6_192_23_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A31 = read_bin('../data/test6_192_31_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A32 = read_bin('../data/test6_192_32_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
A33 = read_bin('../data/test6_192_33_h.dat', dim1 / 3 + 1, dim1 / 3 + 1);
xa_track = zeros(1, length(A11(1,1,:)));
ya_track = zeros(1, length(A11(1,1,:)));

B = read_bin('../data/test6_192_h.dat', dim1 + 1, dim1 + 1);
C = read_bin('../data/test6_384_curl.dat', dim2 + 1, dim2 + 1);
xb_track = zeros(1, length(B(1,1,:)));
yb_track = zeros(1, length(B(1,1,:)));
xc_track = zeros(1, length(C(1,1,:)));
yc_track = zeros(1, length(C(1,1,:)));

xB = [0 : 1 : dim1];
yB = [0 : 1 : dim1];
[XB, YB] = meshgrid(xB, yB);

xC = [0 : dim1 / dim2 : dim1];
yC = [0 : dim1 / dim2 : dim1];
[XC, YC] = meshgrid(xC, yC);


axis([0 dim1 0 dim1 0 10000])
x1 = [0 : 1 : dim1 / 3];
x2 = [dim1 / 3 : 1 : 2 * dim1 / 3];
x3 = [2 * dim1 / 3 : 1 : dim1];
y1 = [0 : 1 : dim1 / 3];
y2 = [dim1 / 3 : 1 : 2 * dim1 / 3];
y3 = [2 * dim1 / 3 : 1 : dim1];

x11 = [0 : dim1 / dim2 : dim1 / 3];
y11 = [0 : dim1 / dim2 : dim1 / 3];
x22 = [dim1 / 3 : dim1 / dim2 : 2 * dim1 / 3];
y22 = [dim1 / 3 : dim1 / dim2 : 2 * dim1 / 3];
x33 = [2 * dim1 / 3 : dim1 / dim2 : dim1];
y33 = [2 * dim1 / 3 : dim1 / dim2 : dim1];

[X11, Y11] = meshgrid(x1, y1);
[X12, Y12] = meshgrid(x1, y2);
[X13, Y13] = meshgrid(x1, y3);
[X21, Y21] = meshgrid(x2, y1);
[X22, Y22] = meshgrid(x22, y22);
[X23, Y23] = meshgrid(x2, y3);
[X31, Y31] = meshgrid(x3, y1);
[X32, Y32] = meshgrid(x3, y2);
[X33, Y33] = meshgrid(x3, y3);

%минимум уровня
hold on; grid on
for i = 1:1:length(C(1,1,:))
    c = C(:,:,i)';
    [m, ind] = min(c, [], 'all');
    plot(XC(ind), YC(ind), '.g')
    c(ind) = 10000;
    [m, ind] = min(c, [], 'all');
    plot(XC(ind), YC(ind), '.g')
end

for i = 1:1:length(A11(1,1,:))
    a11 = A11(:,:,i)';
    a12 = A12(:,:,i)';
    a13 = A13(:,:,i)';
    a21 = A21(:,:,i)';
    a22 = A22(:,:,i)';
    a23 = A23(:,:,i)';
    a31 = A31(:,:,i)';
    a32 = A32(:,:,i)';
    a33 = A33(:,:,i)';

    [mina11, I11] = min(a11, [], 'all');
    [mina12, I12] = min(a12, [], 'all');
    [mina13, I13] = min(a13, [], 'all');
    [mina21, I21] = min(a21, [], 'all');
    [mina22, I22] = min(a22, [], 'all');
    [mina23, I23] = min(a23, [], 'all');
    [mina31, I31] = min(a31, [], 'all');
    [mina32, I32] = min(a32, [], 'all');
    [mina33, I33] = min(a33, [], 'all');

    [m, min_cells] = min([mina11, mina12, mina13, mina21, mina22, mina23, mina31, mina32, mina33]);

    if(min_cells == 1)
        plot(X11(I11), Y11(I11), '.b')
        a11(I11) = 10000;
    elseif (min_cells == 2)
        plot(X12(I12), Y12(I12), '.b')
        a12(I12) = 10000;
    elseif (min_cells == 3)
        plot(X13(I13), Y13(I13), '.b')
        a13(I13) = 10000;
    elseif (min_cells == 4)
        plot(X21(I21), Y21(I21), '.b')
        a21(I21) = 10000;
    elseif (min_cells == 5)
        plot(X22(I22), Y22(I22), '.b')
        a22(I22) = 10000;
    elseif (min_cells == 6)
        plot(X23(I23), Y23(I23), '.b')
        a23(I23) = 10000;
    elseif (min_cells == 7)
        plot(X31(I31), Y31(I31), '.b')
        a31(I31) = 10000;
    elseif (min_cells == 8)
        plot(X32(I32), Y32(I32), '.b')
        a32(I32) = 10000;
    elseif (min_cells == 9)
        plot(X33(I33), Y33(I33), '.b')
        a33(I33) = 10000;
    end

    [mina11, I11] = min(a11, [], 'all');
    [mina12, I12] = min(a12, [], 'all');
    [mina13, I13] = min(a13, [], 'all');
    [mina21, I21] = min(a21, [], 'all');
    [mina22, I22] = min(a22, [], 'all');
    [mina23, I23] = min(a23, [], 'all');
    [mina31, I31] = min(a31, [], 'all');
    [mina32, I32] = min(a32, [], 'all');
    [mina33, I33] = min(a33, [], 'all');

    [m, min_cells] = min([mina11, mina12, mina13, mina21, mina22, mina23, mina31, mina32, mina33]);

    if(min_cells == 1)
        plot(X11(I11), Y11(I11), '.b')
        a11(I11) = 10000;
    elseif (min_cells == 2)
        plot(X12(I12), Y12(I12), '.b')
        a12(I12) = 10000;
    elseif (min_cells == 3)
        plot(X13(I13), Y13(I13), '.b')
        a13(I13) = 10000;
    elseif (min_cells == 4)
        plot(X21(I21), Y21(I21), '.b')
        a21(I21) = 10000;
    elseif (min_cells == 5)
        plot(X22(I22), Y22(I22), '.b')
        a22(I22) = 10000;
    elseif (min_cells == 6)
        plot(X23(I23), Y23(I23), '.b')
        a23(I23) = 10000;
    elseif (min_cells == 7)
        plot(X31(I31), Y31(I31), '.b')
        a31(I31) = 10000;
    elseif (min_cells == 8)
        plot(X32(I32), Y32(I32), '.b')
        a32(I32) = 10000;
    elseif (min_cells == 9)
        plot(X33(I33), Y33(I33), '.b')
        a33(I33) = 10000;
    end
end

for i = 1:1:length(B(1,1,:))
    b = B(:,:,i)';
    [m, ind] = min(b, [], 'all');
    plot(XB(ind), YB(ind), '.-r')
    b(ind) = 10000;
    [m, ind] = min(b, [], 'all');
    plot(XB(ind), YB(ind), '.r')
end
title('MIN')



%взвешенный центр вихревого ядра
figure;
hold on; grid on;
Up_int   = zeros(1, length(A11(1,1,:)));
Down_int = zeros(1, length(A11(1,1,:)));

for t = 1:1:length(A11(1,1,:))
    Down_int(t) = trapz(y3, trapz(x3, A33(:, :, t)));
    Down_int(t) = Down_int(t) + trapz(y2, trapz(x3, A32(:, :, t)));
    Down_int(t) = Down_int(t) + trapz(y1, trapz(x3, A31(:, :, t)));

    Down_int(t) = Down_int(t) + trapz(y3, trapz(x2((dim1 / 3) / 2 + 1 : end : end), A23((dim1 / 3) / 2 + 1 : end, :, t)));
    Down_int(t) = Down_int(t) + trapz(y22, trapz(x22((dim2 / 3) / 2 + 1 : end), A22((dim2 / 3) / 2 + 1 : end, :, t)));
    Down_int(t) = Down_int(t) + trapz(y1, trapz(x2((dim1 / 3) / 2 + 1 : end), A21((dim1 / 3) / 2 + 1 : end, :, t)));


end