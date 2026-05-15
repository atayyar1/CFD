function T = solveGaussSeidel(aE,aW,aN,aS,aC,bC,Nx,Ny,T,maxIter,tol)

for iter = 1:maxIter

    Told = T;

    for i = 1:Nx
    for j = 1:Ny

        TE = 0; TW = 0; TN = 0; TS = 0;

        if i < Nx
            TE = T(i+1,j);
        end
        if i > 1
            TW = T(i-1,j);
        end

        if j < Ny
            TN = T(i,j+1);
        end

        if j > 1
            TS = T(i,j-1);
        end

        T(i,j) = ( bC(i,j) - aE(i,j)*TE - aW(i,j)*TW - aN(i,j)*TN - aS(i,j)*TS ) / aC(i,j);

    end
    end

    % convergence check
    err = max(abs(T(:)-Told(:)));

    if err < tol
        fprintf('Gauss-Seidel converged in %d iterations\n',iter)
        break
    end

end

end
