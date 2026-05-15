clc
clear
close all 
rho = 0.8;          % density          [kg/m^3]
mu  = 5e-5;         % dynamic viscosity [Pa.s]
cp  = 1030;         % specific heat    [J/kg.K]
k   = 0.036;        % thermal cond.    [W/m.K]
  Lx = 2.0;           % domain length [m]
Ly = 1.0;           % domain height [m]
tol_cont = 1e-8;
tol_vel  = 1e-4;

tol_bDC  = 1e-2;
%% --- Boundary condition values ---
T_top    = 300;     % top wall temperature         [K]
T_bot    = 500;     % bottom wall temperature      [K]
T_inlet    = 350;     % inlet temperature            [K]  (professor's note)
v_inlet  = -5.0;    % inlet velocity (downward)    [m/s]
h_left      = 15;      % Robin h on upper-left wall   [W/m^2.K]
T_inf      = 300;     % Robin ambient temperature    [K]
P_ref    = 1e5;     % reference pressure at outlet [Pa]
 
%% --- Inlet / outlet geometry ---
x_inlet_start = 1.8;   % inlet starts 0.2 m from right edge [m]
y_outlet_end  = 0.1;   % outlet occupies bottom 0.1 m of left wall [m]

%% --- SIMPLE solver settings ---
alpha_u   = 0.7;        % velocity under-relaxation     (u)
alpha_v   = 0.7;        % velocity under-relaxation (v)
alpha_p   = 0.3;        % pressure under-relaxation
tol_SIMPLE = 1e-8;    % was 1e-5
maxIter    = 10000;    % keep this, will need more iterations now

 solver.type    = 'TDMA';
solver.maxIter = 5000;
solver.tol     = 1e-8;

 IntFactor    = "volume";
Decomp_method = "minimum";

 %  GRID LOOP
 grid_configs = [80 40];
results = struct();

for g = 1:size(grid_configs,1)

    Nx = grid_configs(g,1);
    Ny = grid_configs(g,2);
    m  = Nx + 1;
    n  = Ny + 1;
     if Nx == 160
    beta_DC = 0.05;
else
    beta_DC = 0.2;
end
     fprintf('  Running grid %d x %d\n', Nx, Ny)
 
    %% --- Grid and geometry ---
    [X,Y]   = GridRect(m, n, Lx, Ly);
    [XC,YC] = getCentroids(X,Y);
    Vc      = getCellVolumes(X,Y);
    [Se,Sw,Sn,Ss]                     = getSurfaceVectors(X,Y,Nx,Ny);
    [Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys]        = getFaceCenters(X,Y,Nx,Ny);
    [CE,CW,CN,CS]                     = getCFVectors(XC,YC,Nx,Ny);
    [ge,gw,gn,gs]                     = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,IntFactor);
    [Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv]    = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,Decomp_method);

     dx = Lx / Nx;
    dy = Ly / Ny;
       %% --- Identify inlet and outlet cell indices ---
    % Inlet: top boundary (j = Ny), cells where x >= x_inlet_start
    i_inlet = find(XC(:,1) > x_inlet_start);

    % Outlet: left boundary (i = 1), cells where y <= y_outlet_end
       j_outlet = find(YC(1,:) < y_outlet_end - dy/4);    % row indices (strict)

    fprintf('  Inlet  cells on top  : i = %d to %d  (%d cells)\n', ...
        i_inlet(1), i_inlet(end), numel(i_inlet))
    fprintf('  Outlet cells on left : j = %d to %d  (%d cells)\n', ...
           j_outlet(1), j_outlet(end), numel(j_outlet))

    %% --- Initialise flow fields ---
    u = zeros(Nx,Ny);       % x-velocity  [m/s]
    v = zeros(Nx,Ny);       % y-velocity  [m/s]
    p = P_ref * ones(Nx,Ny); % pressure   [Pa]

    % Face mass fluxes (initialised to zero)
    Fe = zeros(Nx,Ny);
    Fw = zeros(Nx,Ny);
    Fn = zeros(Nx,Ny);
    Fs = zeros(Nx,Ny);
    bDC_u = zeros(Nx, Ny);
    bDC_v = zeros(Nx, Ny);
    fprintf('  Starting SIMPLE iterations\n')
     res_cont_hist = zeros(maxIter, 1);
    res_u_hist    = zeros(maxIter, 1);
    res_v_hist    = zeros(maxIter, 1);
    res_cont_hist = zeros(maxIter, 1);
    res_u_hist    = zeros(maxIter, 1);
    res_v_hist    = zeros(maxIter, 1);
    smart_active  = false;   
    phase1_done   = false;

% --- pre-allocate convergence history (trimmed after loop) ---
    res_cont_hist = zeros(maxIter, 1);
    res_u_hist    = zeros(maxIter, 1);
    res_v_hist    = zeros(maxIter, 1);

    % --- two-phase SMART flags ---
    smart_active      = false;
    phase1_done       = false;
    phase2_start_iter = 0;

for iter = 1:maxIter

    u_prev    = u;
    v_prev    = v;
    bDC_u_prev = bDC_u;
    bDC_v_prev = bDC_v;

    % ------------------------------------------------
    % 1.  U-MOMENTUM
    % ------------------------------------------------
    [aE_u, aW_u, aN_u, aS_u, aC_u, bC_u, ~] = ...
        assembleMomentumCoeffs('u', u, rho, mu, p, ...
            Fe, Fw, Fn, Fs, ...
            Ee, Ew, En, Es, CE, CW, CN, CS, ...
            Se, Sw, Sn, Ss, Vc, Nx, Ny, alpha_u);

    [aE_u, aW_u, aN_u, aS_u, aC_u, bC_u] = ...
        applyMomentumBC('u', aE_u, aW_u, aN_u, aS_u, aC_u, bC_u, ...
            mu, rho, p, P_ref, ...
            XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, ...
            Se, Sw, Sn, Ss, Fe, Fw, Fn, Fs, ...
            Nx, Ny, i_inlet, j_outlet, v_inlet, dx, dy);

    aP_u = aC_u;

    bC_u = bC_u + (1 - alpha_u) / alpha_u * aP_u .* u;
    aC_u = aC_u / alpha_u;

    if smart_active
        bC_u = bC_u + bDC_u;
    end

    u_star = solveLinearSystem(aE_u, aW_u, aN_u, aS_u, aC_u, bC_u, Nx, Ny, u, solver);

    % ------------------------------------------------
    % 2.  V-MOMENTUM
    % ------------------------------------------------
    [aE_v, aW_v, aN_v, aS_v, aC_v, bC_v, ~] = ...
        assembleMomentumCoeffs('v', v, rho, mu, p, ...
            Fe, Fw, Fn, Fs, ...
            Ee, Ew, En, Es, CE, CW, CN, CS, ...
            Se, Sw, Sn, Ss, Vc, Nx, Ny, alpha_v);

    [aE_v, aW_v, aN_v, aS_v, aC_v, bC_v] = ...
        applyMomentumBC('v', aE_v, aW_v, aN_v, aS_v, aC_v, bC_v, ...
            mu, rho, p, P_ref, ...
            XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, ...
            Se, Sw, Sn, Ss, Fe, Fw, Fn, Fs, ...
            Nx, Ny, i_inlet, j_outlet, v_inlet, dx, dy);

    aP_v = aC_v;

    bC_v = bC_v + (1 - alpha_v) / alpha_v * aP_v .* v;
    aC_v = aC_v / alpha_v;

    if smart_active
        bC_v = bC_v + bDC_v;
    end

    v_star = solveLinearSystem(aE_v, aW_v, aN_v, aS_v, aC_v, bC_v, Nx, Ny, v, solver);

    % ------------------------------------------------
    % 3.  RHIE-CHOW
    % ------------------------------------------------
    [Fe, Fw, Fn, Fs] = computeRhieChow( ...
        u_star, v_star, u, v, ...
        aP_u, aP_v, p, rho, Vc, ...
        ge, gn, Se, Sw, Sn, Ss, ...
        Fe, Fn, ...
        Nx, Ny, dx, dy, ...
        i_inlet, j_outlet, v_inlet, P_ref, ...
        alpha_u, alpha_v);

    % ------------------------------------------------
    % 4.  PRESSURE CORRECTION
    % ------------------------------------------------
    [aE_p, aW_p, aN_p, aS_p, aC_p, bC_p] = ...
        assemblePressureCorrection( ...
            Fe, Fw, Fn, Fs, aP_u, aP_v, Vc, rho, ...
            ge, gn, Nx, Ny, dx, dy, j_outlet);

    p_prime = solveLinearSystem(aE_p, aW_p, aN_p, aS_p, aC_p, bC_p, ...
        Nx, Ny, zeros(Nx,Ny), solver);

    % ------------------------------------------------
    % 5.  CORRECT
    % ------------------------------------------------
    [u, v, p, Fe, Fw, Fn, Fs] = correctSIMPLE( ...
        u_star, v_star, p, p_prime, ...
        aP_u, aP_v, Vc, rho, ...
        Fe, Fw, Fn, Fs, ...
        ge, gn, Nx, Ny, dx, dy, ...
        alpha_p, j_outlet);

    % ------------------------------------------------
    % 6.  UPDATE SMART CORRECTION
    % ------------------------------------------------
    Fe_s = Fe;  Fe_s(Nx, :) = 0;
    Fw_s = Fw;  Fw_s(1,  :) = 0;
    Fn_s = Fn;  Fn_s(:, Ny) = 0;
    Fs_s = Fs;  Fs_s(:,  1) = 0;

    % SMART for u
    [ue_u, uw_u, un_u, us_u] = interpolateTemperatureToFaces( ...
        u, ge, gw, gn, gs, Nx, Ny);
    for i = 1:Nx
        un_u(i, Ny) = 0;
        us_u(i, 1)  = 0;
    end
    for j = 1:Ny
        ue_u(Nx, j) = 0;
        if ~any(j_outlet == j)
            uw_u(1, j) = 0;
        end
    end
    [dudx, dudy] = computeCellGradient( ...
        ue_u, uw_u, un_u, us_u, Se, Sw, Sn, Ss, Vc, Nx, Ny);
    [dudx, dudy] = limitGradients(u, dudx, dudy, XC, YC, Nx, Ny);
    bDC_u_new = computeSMARTCorrection(u, Fe_s, Fw_s, Fn_s, Fs_s, dudx, dudy, ...
        XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, Nx, Ny);

    % SMART for v
    [ve_v, vw_v, vn_v, vs_v] = interpolateTemperatureToFaces( ...
        v, ge, gw, gn, gs, Nx, Ny);
    for i = 1:Nx
        if any(i_inlet == i)
            vn_v(i, Ny) = v_inlet;
        else
            vn_v(i, Ny) = 0;
        end
        vs_v(i, 1) = 0;
    end
    for j = 1:Ny
        ve_v(Nx, j) = 0;
        if ~any(j_outlet == j)
            vw_v(1, j) = 0;
        end
    end
    [dvdx, dvdy] = computeCellGradient( ...
        ve_v, vw_v, vn_v, vs_v, Se, Sw, Sn, Ss, Vc, Nx, Ny);
    [dvdx, dvdy] = limitGradients(v, dvdx, dvdy, XC, YC, Nx, Ny);
    bDC_v_new = computeSMARTCorrection(v, Fe_s, Fw_s, Fn_s, Fs_s, dvdx, dvdy, ...
        XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, Nx, Ny);

    % Blended update
    bDC_u = (1 - beta_DC) * bDC_u + beta_DC * bDC_u_new;
    bDC_v = (1 - beta_DC) * bDC_v + beta_DC * bDC_v_new;

    % ------------------------------------------------
    % 7.  Residuals
    % ------------------------------------------------
    res_cont = sum(abs(Fe + Fw + Fn + Fs), 'all') / (Nx*Ny);
    res_u    = max(abs(u(:) - u_prev(:))) / (max(abs(u(:))) + eps);
    res_v    = max(abs(v(:) - v_prev(:))) / (max(abs(v(:))) + eps);

    % bDC convergence: how much did the SMART correction change this iter?
    bDC_mag  = max(max(abs(bDC_u(:))), max(abs(bDC_v(:)))) + eps;
    bDC_chg  = max(max(abs(bDC_u(:) - bDC_u_prev(:))), ...
                   max(abs(bDC_v(:) - bDC_v_prev(:))));
    res_bDC  = bDC_chg / bDC_mag;

    res_cont_hist(iter) = res_cont;
    res_u_hist(iter)    = res_u;
    res_v_hist(iter)    = res_v;

    if mod(iter, 100) == 0 || iter == 1
        if smart_active
            fprintf('    iter %4d [SMART] | res_cont = %.3e | res_u = %.3e | res_v = %.3e | res_bDC = %.3e\n', ...
                iter, res_cont, res_u, res_v, res_bDC)
        else
            fprintf('    iter %4d [upwind] | res_cont = %.3e\n', iter, res_cont)
        end
    end

    % ------------------------------------------------
    % 8.  Phase transitions and convergence
    % ------------------------------------------------

    % Phase 1 → Phase 2 transition
    if ~phase1_done && res_cont < tol_SIMPLE
        phase1_done       = true;
        smart_active      = true;
        phase2_start_iter = iter + 1;
        fprintf('  >>> Phase 1 (upwind) converged at iter %d  (res_cont = %.3e)\n', ...
            iter, res_cont)
        fprintf('  >>> Activating SMART — starting Phase 2\n')
        bDC_u = zeros(Nx, Ny);
        bDC_v = zeros(Nx, Ny);
        continue
    end

    % Phase 2 convergence: continuity tight AND SMART correction settled
    if smart_active && res_cont < tol_SIMPLE && res_bDC < tol_bDC
        fprintf('  Converged (SMART) at iter %d  (res_cont = %.3e)\n', iter, res_cont)
        fprintf('    res_u = %.3e | res_v = %.3e | res_bDC = %.3e\n', ...
            res_u, res_v, res_bDC)
        break
    end

    if iter == maxIter
        fprintf('  WARNING: did not converge in %d iterations\n', maxIter)
        fprintf('    smart_active = %d\n', smart_active)
        fprintf('    res_cont = %.3e | res_u = %.3e | res_v = %.3e | res_bDC = %.3e\n', ...
            res_cont, res_u, res_v, res_bDC)
    end

end   % SIMPLE loop

    % Trim history
    res_cont_hist = res_cont_hist(1:iter);
    res_u_hist    = res_u_hist(1:iter);
    res_v_hist    = res_v_hist(1:iter);

     res_cont_hist = res_cont_hist(1:iter);
    res_u_hist    = res_u_hist(1:iter);
    res_v_hist    = res_v_hist(1:iter);

     res_cont_hist = res_cont_hist(1:iter);
    res_u_hist    = res_u_hist(1:iter);
    res_v_hist    = res_v_hist(1:iter);

     fprintf('  Solving energy equation...\n')
    T = solveEnergyA4( ...
        u, v, Fe, Fw, Fn, Fs, ...
        k, rho, cp, ...
        Ee, Ew, En, Es, CE, CW, CN, CS, ...
        Tev, Twv, Tnv, Tsv, ...
        ge, gw, gn, gs, ...
        Se, Sw, Sn, Ss, Vc, ...
        XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, ...
        Nx, Ny, ...
        i_inlet, j_outlet, ...
        T_top, T_bot, T_inlet, h_left, T_inf, ...
        solver);

    %% --- Store results ---
    results(g).Nx           = Nx;
    results(g).Ny           = Ny;
    results(g).XC           = XC;
    results(g).YC           = YC;
    results(g).u            = u;
    results(g).v            = v;
    results(g).p            = p;
    results(g).T            = T;
    results(g).res_cont_hist = res_cont_hist;
    results(g).res_u_hist    = res_u_hist;
    results(g).res_v_hist    = res_v_hist;

    fprintf('  Grid %dx%d complete.\n', Nx, Ny)

end    
for g = 1:size(grid_configs,1)
figure;
semilogy(results(g).res_cont_hist, 'b-',  'LineWidth', 1.5); hold on
semilogy(results(g).res_u_hist,    'r--', 'LineWidth', 1.5)
semilogy(results(g).res_v_hist,    'g--', 'LineWidth', 1.5)
yline(tol_SIMPLE, 'k:', 'LineWidth', 1)
xlabel('Iteration'); ylabel('Residual')
legend('Continuity', 'u-velocity', 'v-velocity', 'Tolerance')
title(sprintf('SIMPLE Convergence — %d\\times%d', results(g).Nx, results(g).Ny))
grid on
    plotFlowResults(results(g).XC, results(g).YC, ...
                    results(g).u,  results(g).v, ...
                    results(g).p,  results(g).T, ...
                    results(g).Nx, results(g).Ny);

end