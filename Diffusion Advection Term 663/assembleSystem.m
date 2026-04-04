function [A,b] = assembleSystem(aE,aW,aN,aS,aC,bC,Nx,Ny)

N = Nx*Ny;

A = zeros(N,N);
b = zeros(N,1);

for i = 1:Nx
    for j = 1:Ny

        % global index
        k = (j-1)*Nx + i;

        % central coefficient
        A(k,k) = aC(i,j);

        % east neighbor
        if i < Nx
            kE = k + 1;
            A(k,kE) = aE(i,j);
        end

        % west neighbor
        if i > 1
            kW = k - 1;
            A(k,kW) = aW(i,j);
        end

        % north neighbor
        if j < Ny
            kN = k + Nx;
            A(k,kN) = aN(i,j);
        end

        % south neighbor
        if j > 1
            kS = k - Nx;
            A(k,kS) = aS(i,j);
        end

        % RHS
        b(k) = bC(i,j);

    end
end

end