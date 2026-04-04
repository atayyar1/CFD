function T = solveTDMA(aE,aW,aN,aS,aC,bC,Nx,Ny,T,maxIter,tol)
for iter = 1:maxIter
    Told = T;
    for j = 1:Ny
        lower = zeros(Nx,1);
        diagMain = zeros(Nx,1);
        upper = zeros(Nx,1);
        rhs = zeros(Nx,1);
        for i = 1:Nx
            TS = 0.0;
            TN = 0.0;
            if j > 1
                TS = T(i,j-1);
            end
            if j < Ny
                TN = T(i,j+1);
            end
            diagMain(i) = aC(i,j);
            rhs(i) = bC(i,j) - aS(i,j)*TS - aN(i,j)*TN;
            if i > 1
                lower(i) = aW(i,j);
            end
            if i < Nx
                upper(i) = aE(i,j);
            end
        end
        T(:,j) = solveTridiagonal(lower,diagMain,upper,rhs);
    end
    err = max(abs(T(:) - Told(:)));
    if err < tol
        fprintf('TDMA converged in %d iterations\n',iter)
        break
    end

end

end

function x = solveTridiagonal(a,b,c,d)
n = length(b);

cp = zeros(n,1);
dp = zeros(n,1);
x = zeros(n,1);

cp(1) = c(1) / b(1);
dp(1) = d(1) / b(1);

for i = 2:n
    denom = b(i) - a(i)*cp(i-1);

    if i < n
        cp(i) = c(i) / denom;
    end

    dp(i) = (d(i) - a(i)*dp(i-1)) / denom;
end

x(n) = dp(n);

for i = n-1:-1:1
    x(i) = dp(i) - cp(i)*x(i+1);
end

end
%% Please note that I used chat gpt to debug this file