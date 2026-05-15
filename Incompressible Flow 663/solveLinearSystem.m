function T = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T,solver)
%SOLVELINEARSYSTEM choices:
%   solver.type = 'gauss-seidel' or 'tdma'

switch lower(solver.type)
    case 'gauss-seidel'
        T = solveGaussSeidel( ...
            aE,aW,aN,aS,aC,bC, ...
            Nx,Ny,T, ...
            solver.maxIter,solver.tol);

    case 'tdma'
        T = solveTDMA( ...
            aE,aW,aN,aS,aC,bC, ...
            Nx,Ny,T, ...
            solver.maxIter,solver.tol);

    otherwise
        error('Unknown solver type. Use gauss-seidel or tdma.');
end

end
