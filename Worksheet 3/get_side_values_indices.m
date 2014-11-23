function [ i, j ] = get_side_values_indices( Nx, Ny )

N = Nx-1;
i = zeros(1, N*N);
j = zeros(1, N*N);

for k=0:N    
    i(1+k*N:(k+1)*N) = ((k+1)+k*N): (k + (k+1)*N);
    j(1+k*N:(k+1)*N) = ((k+2)+k*N): (k + 1 + (k+1)*N);
end

end

