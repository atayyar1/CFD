function T_anal = computeAnalyticalSolution(XC, YC, t, a, b, k, rho, cp)

% because of the Neumann BCs in x, the x eigenfunctions are cos(m*pi*x/a)
% but since IC is uniform in x, all m>=1 terms vanish (integral of cos over
% full period = 0), so solution only depends on y

alpha = k / (rho * cp);

[Nx, Ny] = size(XC);
T_anal = 1000*zeros(Nx, Ny);
N_terms = 100;
for n = 1:N_terms
    An = (2/b) * 1000 * (b/(n*pi)) * (1 - cos(n*pi));    
    lambda_n = (n*pi/b)^2;
    T_anal = T_anal + An * sin(n*pi*YC/b) * exp(-lambda_n * alpha * t);
end
end