function [A] = read_bin(file, Nx, Ny)
    id = fopen(file);
    A=fread(id, 'single');
    A=reshape(A, Nx, Ny, []);
end

