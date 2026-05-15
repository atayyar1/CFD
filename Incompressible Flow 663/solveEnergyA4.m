function T = solveEnergyA4(u, v, Fe, Fw, Fn, Fs, k, rho, cp, ...
    Ee, Ew, En, Es, CE, CW, CN, CS,  Tev, Twv, Tnv, Tsv, ge, gw, gn, gs, Se, Sw, Sn, Ss, Vc, ...
    XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, ...
    Nx, Ny, i_inlet, j_outlet,T_top, T_bot, T_inlet, h_left, T_inf, solver)

Gamma_e = k* ones(Nx, Ny);
Gamma_w = k*ones(Nx, Ny);
Gamma_n = k*ones(Nx, Ny);
Gamma_s = k*ones(Nx, Ny);
S_zero= zeros(Nx, Ny);
T = (T_top + T_bot)/2*ones(Nx, Ny);
maxDC = 50;
tolDC = 1e-8;
for dc = 1:maxDC
    T_old = T;
    [Te, Tw, Tn, Ts]=interpolateTemperatureToFaces(T, ge, gw, gn, gs, Nx, Ny); 
    for i = 1:Nx
        if any(i_inlet == i)
            Tn(i,Ny) = T_inlet;
        else
        Tn(i,Ny) = T_top;
        end
    end
    for i = 1:Nx
    Ts(i,1) = T_bot;
    end
     for j = 1:Ny
        Te(Nx,j) = T(Nx,j);
    end
     for j = 1:Ny
        if any(j_outlet == j)
             Tw(1,j) = T(1,j);
        else
             dist_w  = abs(XC(1,j) - Xw(1,j));
            D_robin = 1 / (dist_w/k + 1/h_left);
             Tw(1,j) = (k/dist_w * T(1,j) + h_left * T_inf) / (k/dist_w + h_left);
        end
    end
    [dTdx, dTdy] = computeCellGradient(Te, Tw, Tn, Ts, Se, Sw, Sn, Ss, Vc, Nx, Ny);
 [dTdx, dTdy] = limitGradients(T, dTdx, dTdy, XC, YC, Nx, Ny);   

    Fe_h = cp * Fe;
    Fw_h = cp * Fw;
    Fn_h = cp * Fn;
    Fs_h = cp * Fs;

     bCorr = computeNonOrthogonalCorrection(Gamma_e, Gamma_w, Gamma_n, Gamma_s, ...
        dTdx, dTdy, Tev, Twv, Tnv, Tsv, ge, gw, gn, gs, Nx, Ny);

     [aE, aW, aN, aS, aC, bC] = computeDiffusionCoefficients( ...
        Gamma_e, Gamma_w, Gamma_n, Gamma_s, ...
        Ee, Ew, En, Es, CE, CW, CN, CS, ...
        S_zero, Vc, bCorr, Nx, Ny);

     [aE, aW, aN, aS, aC] = assembleConvectionUpwind( ...
        aE, aW, aN, aS, aC, Fe_h, Fw_h, Fn_h, Fs_h, Nx, Ny);

 
   for i = 1:Nx
        area_n = norm(squeeze(Sn(i,Ny,:)));
        dist_n = abs(Yn(i,Ny) - YC(i,Ny));
        D_n    = k * area_n / dist_n;

        if any(i_inlet == i)
            T_b = T_inlet;   
        else
            T_b = T_top;        
        end

        F_b = Fn_h(i,Ny);    
        aN(i,Ny) = 0;
        aC(i,Ny) = aC(i,Ny) + D_n + max(F_b, 0);
        bC(i,Ny) = bC(i,Ny) + (D_n + max(-F_b, 0)) * T_b;
    end

     for i = 1:Nx
        area_s = norm(squeeze(Ss(i,1,:)));
        dist_s = abs(YC(i,1) - Ys(i,1));
        D_s    = k * area_s / dist_s;
        F_b    = 0;  

        aS(i,1) = 0;
        aC(i,1) = aC(i,1) + D_s;
        bC(i,1) = bC(i,1) + D_s * T_bot;
    end
   for j = 1:Ny
        aE(Nx,j) = 0;
     end

     for j = 1:Ny
        area_w = norm(squeeze(Sw(1,j,:)));
        dist_w = abs(XC(1,j) - Xw(1,j));
        aW(1,j) = 0;

        if any(j_outlet == j)
        else
             D_robin    = area_w / (dist_w/k + 1/h_left);
            F_b       = 0;    
            aC(1,j) = aC(1,j) + D_robin;
            bC(1,j) = bC(1,j) + D_robin * T_inf;
        end
    end
    bDC = computeSMARTCorrection(T, Fe_h, Fw_h, Fn_h, Fs_h, dTdx, dTdy, ...
        XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, Nx, Ny);
    bC = bC + bDC;
    T = solveLinearSystem(aE, aW, aN, aS, aC, bC, Nx, Ny, T, solver);

     if max(abs(T(:) - T_old(:))) < tolDC
        fprintf('    Energy: converged in %d deferred correction iterations\n', dc)
        break
    end
end

end