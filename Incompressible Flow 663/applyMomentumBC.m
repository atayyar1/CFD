function [aE, aW, aN, aS, aC, bC] = applyMomentumBC( comp, aE, aW, aN, aS, aC, bC, mu, rho, p, P_ref, ...
    XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, ...
    Se, Sw, Sn, Ss, Fe, Fw, Fn, Fs,  Nx, Ny, i_inlet, j_outlet, v_inlet, dx, dy)

for i = 1:Nx
    area_n = norm(squeeze(Sn(i,Ny,:)));
    dist_n = abs(Yn(i,Ny) - YC(i,Ny));
    D_n    = mu * area_n / dist_n;

    if any(i_inlet == i)
        switch lower(comp)
            case 'u';  phi_b = 0;
            case 'v';  phi_b = v_inlet;
        end
        F_b = rho * v_inlet * area_n;   % known inlet flux, not stale Fn
    else
          phi_b = 0;
        F_b   = 0;
    end

    aN(i,Ny) = 0;
    aC(i,Ny) = aC(i,Ny) + D_n + max(F_b, 0);
    bC(i,Ny) = bC(i,Ny) + (D_n + max(-F_b, 0)) * phi_b;
end

% 2. south
for i = 1:Nx

    area_s = norm(squeeze(Ss(i,1,:)));
    dist_s = abs(YC(i,1) - Ys(i,1));
    D_s    = mu * area_s / dist_s;

    % Zero south neighbour (may have wrong upwind term)
    aS(i,1) = 0;

    % phi_b = 0 for no-slip → bC term vanishes; only diffusion added to diagonal
    F_b = 0;   % wall: no flux
    aC(i,1) = aC(i,1) + D_s + max(F_b, 0);
 
end

% 3. east
for j = 1:Ny

    area_e = norm(squeeze(Se(Nx,j,:)));
    dist_e = abs(Xe(Nx,j) - XC(Nx,j));
    D_e    = mu * area_e / dist_e;

    aE(Nx,j) = 0;

    F_b = 0;   % wall
    aC(Nx,j) = aC(Nx,j) + D_e + max(F_b, 0);
 
end

% 4. west 
for j = 1:Ny

    area_w = norm(squeeze(Sw(1,j,:)));
    dist_w = abs(XC(1,j) - Xw(1,j));
    D_w    = mu * area_w / dist_w;

    if any(j_outlet == j)
        
        aW(1,j) = 0;   % explicit (should already be 0)

        if strcmpi(comp, 'u') 
             bC(1,j) = bC(1,j) + (P_ref - p(1,j)) * area_w;
        end
        % v-momentum: no correction (west face surface vector has zero y-component)

    else
         aW(1,j) = 0;   

        F_b = 0;
        aC(1,j) = aC(1,j) + D_w + max(F_b, 0);
        % bC(1,j) += 0  (phi_b = 0)

    end
end

end