function [Fe, Fw, Fn, Fs] = computeRhieChow(u_star, v_star, u_old, v_old, ...
    aP_u, aP_v, p, rho, Vc, ...
    ge, gn, Se, Sw, Sn, Ss,Fe_old, Fn_old,Nx, Ny, dx, dy, i_inlet, j_outlet, v_inlet, P_ref, alpha_u, alpha_v)
Fe = zeros(Nx, Ny);
Fw = zeros(Nx, Ny);
Fn = zeros(Nx, Ny);
Fs = zeros(Nx, Ny);

area_e = dy;    
area_n = dx;    

dpdx = zeros(Nx, Ny);
dpdy = zeros(Nx, Ny);

for j = 1:Ny
    for i = 1:Nx
        if i > 1 && i < Nx
            dpdx(i,j) = (p(i+1,j) - p(i-1,j)) / (2*dx);
        elseif i == 1
            dpdx(i,j) = (p(2,j)   - p(1,j))   / dx;
        else  
            dpdx(i,j) = (p(Nx,j)  - p(Nx-1,j)) / dx;
        end

        % y-gradient
        if j > 1 && j < Ny
            dpdy(i,j) = (p(i,j+1) - p(i,j-1)) / (2*dy);
        elseif j == 1
            dpdy(i,j) = (p(i,2)   - p(i,1))   / dy;
        else  % j == Ny
            dpdy(i,j) = (p(i,Ny)  - p(i,Ny-1)) / dy;
        end
    end
end
 
Du = Vc ./ aP_u;   % [Nx,Ny]
Dv = Vc ./ aP_v;
 
for i = 1:Nx-1
    for j = 1:Ny

        gC = ge(i,j);           
        gF = 1 - gC;         

         u_bar = gC*u_star(i,j) + gF*u_star(i+1,j);

         D_bar = gC*Du(i,j) + gF*Du(i+1,j);

         dp_compact = (p(i+1,j) - p(i,j)) / dx;

         dp_interp = gC*dpdx(i,j) + gF*dpdx(i+1,j);
        uf_e_old  = Fe_old(i,j) / (rho * area_e);
        u_bar_old = gC*u_old(i,j) + gF*u_old(i+1,j);
        ur_corr   = (1 - alpha_u) * (uf_e_old - u_bar_old);

        % --- Rhie-Chow face velocity ---
        uf_e = u_bar - D_bar*(dp_compact - dp_interp) + ur_corr;

        % --- face mass flux (outward positive for cell C) ---
        Fe(i,j) = rho * uf_e * area_e;

    end
end
 
for i = 1:Nx
    for j = 1:Ny-1

        gC = gn(i,j);
        gF = 1 - gC;

         v_bar = gC*v_star(i,j) + gF*v_star(i,j+1);

         D_bar = gC*Dv(i,j) + gF*Dv(i,j+1);

         dp_compact = (p(i,j+1) - p(i,j)) / dy;

         dp_interp = gC*dpdy(i,j) + gF*dpdy(i,j+1);

            vf_n_old  = Fn_old(i,j) / (rho * area_n);
        v_bar_old = gC*v_old(i,j) + gF*v_old(i,j+1);
        ur_corr   = (1 - alpha_v) * (vf_n_old - v_bar_old);

        % --- Rhie-Chow face velocity ---
        vf_n = v_bar - D_bar*(dp_compact - dp_interp) + ur_corr;

        % --- face mass flux ---
        Fn(i,j) = rho * vf_n * area_n;

    end
end
for i = 2:Nx
    for j = 1:Ny
        Fw(i,j) = -Fe(i-1,j);
    end
end

for i = 1:Nx
    for j = 2:Ny
        Fs(i,j) = -Fn(i,j-1);
    end
end 
for i = 1:Nx
    if any(i_inlet == i)
        Fn(i,Ny) = rho * v_inlet * area_n; 
    else
        Fn(i,Ny) = 0;
    end
end
 
for i = 1:Nx
    Fs(i,1) = 0;
end
 
for j = 1:Ny
    Fe(Nx,j) = 0;
end
 
for j = 1:Ny
    if any(j_outlet == j)
        dp_boundary = (p(1,j) - P_ref) * 2/dx;   % (∂p/∂x) at boundary face
        u_b = u_star(1,j) - Du(1,j) * (dp_boundary - dpdx(1,j));

        Fw(1,j) = -rho * u_b * area_e;   
    else
        Fw(1,j) = 0;   
    end
end

end