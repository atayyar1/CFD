function [dphidx_lim, dphidy_lim] = limitGradients(phi, dphidx, dphidy, XC, YC, Nx, Ny)
 
dphidx_lim = dphidx;
dphidy_lim = dphidy;

for i = 1:Nx
    for j = 1:Ny
        phi_C = phi(i,j);
         phi_max = phi_C;
        phi_min = phi_C;
        if i < Nx
            phi_max = max(phi_max, phi(i+1,j));
            phi_min = min(phi_min, phi(i+1,j));
        end
        if i > 1
            phi_max = max(phi_max, phi(i-1,j));
            phi_min = min(phi_min, phi(i-1,j));
        end
        if j < Ny
            phi_max = max(phi_max, phi(i,j+1));
            phi_min = min(phi_min, phi(i,j+1));
        end
        if j > 1
            phi_max = max(phi_max, phi(i,j-1));
            phi_min = min(phi_min, phi(i,j-1));
        end
         psi = 1.0;
        neighbors_dx = [];
        neighbors_dy = [];
        if i < Nx
            neighbors_dx(end+1) = XC(i+1,j) - XC(i,j);
            neighbors_dy(end+1) = YC(i+1,j) - YC(i,j);
        end
        if i > 1
            neighbors_dx(end+1) = XC(i-1,j) - XC(i,j);
            neighbors_dy(end+1) = YC(i-1,j) - YC(i,j);
        end
        if j < Ny
            neighbors_dx(end+1) = XC(i,j+1) - XC(i,j);
            neighbors_dy(end+1) = YC(i,j+1) - YC(i,j);
        end
        if j > 1
            neighbors_dx(end+1) = XC(i,j-1) - XC(i,j);
            neighbors_dy(end+1) = YC(i,j-1) - YC(i,j);
        end

        for k = 1:numel(neighbors_dx)
             dphi_pred = dphidx(i,j) * neighbors_dx(k) + ...
                        dphidy(i,j) * neighbors_dy(k);

            if dphi_pred > 1e-14
                 psi_k = min(1.0, (phi_max - phi_C) / dphi_pred);
            elseif dphi_pred < -1e-14
                 psi_k = min(1.0, (phi_min - phi_C) / dphi_pred);
            else
                 psi_k = 1.0;
            end

            psi = min(psi, psi_k);
        end

        psi = max(psi, 0.0);

        dphidx_lim(i,j) = psi * dphidx(i,j);
        dphidy_lim(i,j) = psi * dphidy(i,j);

    end
end

end