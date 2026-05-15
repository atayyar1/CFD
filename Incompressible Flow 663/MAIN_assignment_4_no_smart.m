clc
clear
close all

%% =========================================================
%  MECH 663 — Assignment 4: 2-D Incompressible Flow (SIMPLE)
%  Domain : 2 m (x) × 1 m (y), uniform rectangular grid
%  Fluid  : Air
%  Grids  : 40×20, 80×40, 160×80
% ==========================================================

%% --- Fluid properties ---
rho = 0.8;          % density          [kg/m^3]
mu  = 5e-5;         % dynamic viscosity [Pa.s]
cp  = 1030;         % specific heat    [J/kg.K]
k   = 0.036;        % thermal cond.    [W/m.K]

%% --- Domain ---
Lx = 2.0;           % domain length [m]
Ly = 1.0;           % domain height [m]

%% --- Boundary condition values ---
T_top    = 300;     % top wall temperature         [K]
T_bot    = 500;     % bottom wall temperature      [K]
T_inlet  = 350;     % inlet temperature            [K]  (professor's note)
v_inlet  = -5.0;    % inlet velocity (downward)    [m/s]
h_left   = 15;      % Robin h on upper-left wall   [W/m^2.K]
T_inf    = 300;     % Robin ambient temperature    [K]
P_ref    = 1e5;     % reference pressure at outlet [Pa]

%% --- Inlet / outlet geometry ---
x_inlet_start = 1.8;   % inlet starts 0.2 m from right edge [m]
y_outlet_end  = 0.1;   % outlet occupies bottom 0.1 m of left wall [m]

%% --- SIMPLE solver settings ---
alpha_u   = 0.7;        % velocity under-relaxation     (u)
alpha_v   = 0.7;        % velocity under-relaxation (v)
alpha_p   = 0.3;        % pressure under-relaxation
tol_SIMPLE = 1e-8;    % was 1e-5
maxIter    = 3000;    % keep this, will need more iterations now

%% --- Linear solver settings ---
solver.type    = 'TDMA';
solver.maxIter = 5000;
solver.tol     = 1e-8;

%% --- Geometry decomposition options (same as previous assignments) ---
IntFactor    = "volume";
Decomp_method = "minimum";

%% =========================================================
%  GRID LOOP
% ==========================================================
grid_configs = [40 20; 80 40; 160 80];

for g = 1:size(grid_configs,1)

    Nx = grid_configs(g,1);
    Ny = grid_configs(g,2);
    m  = Nx + 1;
    n  = Ny + 1;

    fprintf('\n========================================\n')
    fprintf('  Running grid %d x %d\n', Nx, Ny)
    fprintf('========================================\n')

    %% --- Grid and geometry ---
    [X,Y]   = GridRect(m, n, Lx, Ly);
    [XC,YC] = getCentroids(X,Y);
    Vc      = getCellVolumes(X,Y);
    [Se,Sw,Sn,Ss]                     = getSurfaceVectors(X,Y,Nx,Ny);
    [Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys]        = getFaceCenters(X,Y,Nx,Ny);
    [CE,CW,CN,CS]                     = getCFVectors(XC,YC,Nx,Ny);
    [ge,gw,gn,gs]                     = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,IntFactor);
    [Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv]    = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,Decomp_method);

    % Cell sizes (uniform grid — used in pressure correction)
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

    % SMART deferred-correction sources for momentum (lagged: updated each
    % SIMPLE iteration using the solution from the previous iteration).
    % Initialised to zero so the first iteration is pure upwind.
    bDC_u = zeros(Nx, Ny);
    bDC_v = zeros(Nx, Ny);

    %% =====================================================
    %  SIMPLE LOOP
    % ======================================================
    fprintf('  Starting SIMPLE iterations...\n')
for iter = 1:maxIter

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

    % Save aP AFTER BCs, BEFORE under-relaxation (correct D for Rhie-Chow)
    aP_u = aC_u;

    % Apply under-relaxation here in MAIN
    bC_u = bC_u + (1 - alpha_u) / alpha_u * aP_u .* u;
    aC_u = aC_u / alpha_u;

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

    % Save aP AFTER BCs, BEFORE under-relaxation
    aP_v = aC_v;

    % Apply under-relaxation
    bC_v = bC_v + (1 - alpha_v) / alpha_v * aP_v .* v;
    aC_v = aC_v / alpha_v; 
    
    v_star = solveLinearSystem(aE_v, aW_v, aN_v, aS_v, aC_v, bC_v, Nx, Ny, v, solver);

    % ------------------------------------------------
    % 3.  RHIE-CHOW face velocities from u*, v*
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
    % 5.  CORRECT velocities, pressure, face fluxes
    % ------------------------------------------------
    [u, v, p, Fe, Fw, Fn, Fs] = correctSIMPLE( ...
        u_star, v_star, p, p_prime, ...
        aP_u, aP_v, Vc, rho, ...
        Fe, Fw, Fn, Fs, ...
        ge, gn, Nx, Ny, dx, dy, ...
        alpha_p, j_outlet);

    % ------------------------------------------------
    % 6.  UPDATE SMART CORRECTIONS FOR NEXT ITERATION
    %     Computed from u, v  (corrected, divergence-free)
    %     and Fe/Fw/Fn/Fs     (corrected, mass-conservative)
    %     NOT from u_star/v_star which are not divergence-free.
    % ------------------------------------------------

    % Mask boundary face fluxes so SMART does not override BCs
    Fe_s = Fe;  Fe_s(Nx, :) = 0;
    Fw_s = Fw;  Fw_s(1,  :) = 0;
    Fn_s = Fn;  Fn_s(:, Ny) = 0;
    Fs_s = Fs;  Fs_s(:,  1) = 0;

    % --- SMART correction for u ---
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
    bDC_u = computeSMARTCorrection(u, Fe_s, Fw_s, Fn_s, Fs_s, dudx, dudy, ...
        XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, Nx, Ny);

    % --- SMART correction for v ---
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
    bDC_v = computeSMARTCorrection(v, Fe_s, Fw_s, Fn_s, Fs_s, dvdx, dvdy, ...
        XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, Nx, Ny);

    % ------------------------------------------------
    % 7.  Convergence check
    % ------------------------------------------------
    res_cont = sum(abs(Fe + Fw + Fn + Fs), 'all') / (Nx*Ny);

    if mod(iter, 100) == 0 || iter == 1
        fprintf('    iter %4d | continuity residual = %.3e\n', iter, res_cont)
    end

    if res_cont < tol_SIMPLE
        fprintf('  Converged at iteration %d  (residual = %.3e)\n', iter, res_cont)
        break
    end

    if iter == maxIter
        fprintf('  WARNING: SIMPLE did not converge in %d iterations\n', maxIter)
    end

end   % SIMPLE loop

    %% =====================================================
    %  ENERGY EQUATION  (solved once with converged u, v)
    % ======================================================
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

    %% =====================================================
    %  STORE RESULTS  (for grid refinement study)
    % ======================================================
    results(g).Nx  = Nx;
    results(g).Ny  = Ny;
    results(g).XC  = XC;
    results(g).YC  = YC;
    results(g).u   = u;
    results(g).v   = v;
    results(g).p   = p;
    results(g).T   = T;

    fprintf('  Grid %dx%d complete.\n', Nx, Ny)

end   % grid loop

%% =========================================================
%  PLOTS  — one figure per grid (velocity, streamlines, T)
%  NOTE: used an LLM to help with the plotting functions
% ==========================================================
for g = 1:size(grid_configs,1)
    plotFlowResults(results(g).XC, results(g).YC, ...
                    results(g).u,  results(g).v, ...
                    results(g).p,  results(g).T, ...
                    results(g).Nx, results(g).Ny);
end