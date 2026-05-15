function [u, v, p, Fe, Fw, Fn, Fs] = correctSIMPLE( ...
    u_star, v_star, p_old, p_prime, ...
    aP_u, aP_v, Vc, rho, ...
    Fe, Fw, Fn, Fs, ...
    ge, gn, Nx, Ny, dx, dy, ...
    alpha_p, j_outlet)
 
Du = Vc ./ aP_u;
Dv = Vc ./ aP_v;

 p = p_old + alpha_p * p_prime;

 
dp_prime_dx = zeros(Nx, Ny);
dp_prime_dy = zeros(Nx, Ny);

for j = 1:Ny
    for i = 1:Nx
        if i > 1 && i < Nx
            dp_prime_dx(i,j) = (p_prime(i+1,j) - p_prime(i-1,j)) / (2*dx);
        elseif i == 1
            dp_prime_dx(i,j) = (p_prime(2,j)  - p_prime(1,j))  / dx;
        else
            dp_prime_dx(i,j) = (p_prime(Nx,j) - p_prime(Nx-1,j)) / dx;
        end

        if j > 1 && j < Ny
            dp_prime_dy(i,j) = (p_prime(i,j+1) - p_prime(i,j-1)) / (2*dy);
        elseif j == 1
            dp_prime_dy(i,j) = (p_prime(i,2)  - p_prime(i,1))  / dy;
        else
            dp_prime_dy(i,j) = (p_prime(i,Ny) - p_prime(i,Ny-1)) / dy;
        end
    end
end
 
u = u_star - Du .* dp_prime_dx;
v = v_star - Dv .* dp_prime_dy;
 
for i = 1:Nx-1
    for j = 1:Ny
        gC  = ge(i,j);
        gF  = 1 - gC;
        D_e = gC*Du(i,j) + gF*Du(i+1,j);

        dm_e      = -rho * D_e * (p_prime(i+1,j) - p_prime(i,j)) / dx * dy;
        Fe(i,j)   = Fe(i,j) + dm_e;
        Fw(i+1,j) = Fw(i+1,j) - dm_e;   % same face, opposite sign for west
    end
end

 for i = 1:Nx
    for j = 1:Ny-1
        gC  = gn(i,j);
        gF  = 1 - gC;
        D_n = gC*Dv(i,j) + gF*Dv(i,j+1);

        dm_n      = -rho * D_n * (p_prime(i,j+1) - p_prime(i,j)) / dy * dx;
        Fn(i,j)   = Fn(i,j) + dm_n;
        Fs(i,j+1) = Fs(i,j+1) - dm_n;   % same face, opposite sign
    end
end
 
for j = 1:Ny
    if any(j_outlet == j)
        dm_outlet = rho * Du(1,j) * p_prime(1,j) * 2/dx * dy;
        Fw(1,j)   = Fw(1,j) + dm_outlet;
    end
end
 
end