clc
clear
close all

%% grid size
m =10; n = 10;
Nx = m-1; Ny = n-1;

%% velocity field - uniform flow, change u and v here if needed
% for 45 degrees: u = 1/sqrt(2), v = 1/sqrt(2)
% for horizontal: u = 1, v = 0 etc
u = 5/sqrt(2);
v = 5/sqrt(2);
rho = 1;

%% solver settings
maxOuterIter = 500;
outerTol = 1e-8;

% inner linear solver
solver.type = 'gauss-seidel';
solver.maxIter = 5000;
solver.tol = 1e-10;

% SMART outer iterations
maxSMARTIter = 100;

%% settings (same as assignment 1)
IntFactor = "volume";
Decomp_method = "minimum";

%% boundary conditions (same as assignment 1)
bc.west.type = 'dirichlet';
bc.west.value = 400.0;

bc.east.type = 'robin';
bc.east.h = 15.0;
bc.east.Tinf = 300.0;

bc.south.type = 'dirichlet';
bc.south.value = 320.0;

bc.north.type = 'neumann';
bc.north.value = 0.0;

%% initial guess
T = 350*ones(Nx,Ny);

%% grid and geometry
[X,Y] = Grid(m,n);
[XC,YC] = getCentroids(X,Y);
Vc = getCellVolumes(X,Y);
[Se,Sw,Sn,Ss] = getSurfaceVectors(X,Y,Nx,Ny);
[Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys] = getFaceCenters(X,Y,Nx,Ny);
[CE,CW,CN,CS] = getCFVectors(XC,YC,Nx,Ny);
[ge,gw,gn,gs] = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,IntFactor);
[Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,Decomp_method);

% face fluxes are constant for uniform flow, compute once
[Fe,Fw,Fn,Fs] = computeFaceFlux(rho,u,v,Se,Sw,Sn,Ss,Nx,Ny);

%% ---- UPWIND SOLVE ----
% nonlinear outer loop (for Gamma and S dependence on T)
upwind_err_hist = [];
for outerIter = 1:maxOuterIter
    Told = T;

    % material properties
    [Gamma,S] = computeMaterialProperties(XC,YC,T,Nx,Ny);

    % gamma and temperature to faces
    [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);
    [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(T,ge,gw,gn,gs,Nx,Ny);
    [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T,Gamma,bc,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);

    % cell gradients and non-orthogonal correction
    [dTdx,dTdy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);
    bCorr =computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);

    % diffusion coefficients
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        Ee,Ew,En,Es,CE,CW,CN,CS,S,Vc,bCorr,Nx,Ny);

    % add convection (upwind)
    [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny);

    % boundary conditions (diffusion BCs, no convective here)
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc);

    % solve
    T = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T,solver);

    outerErr = max(abs(T(:) - Told(:)));
    upwind_err_hist(end+1) = outerErr;
    fprintf('Upwind iter %d, max dT = %e\n', outerIter, outerErr);

    if outerErr < outerTol
        fprintf('Upwind converged in %d iterations\n', outerIter);
        break
    end
end

T_upwind = T;

%% ---- SMART SOLVE ----
% start from upwind solution
T = 350*ones(Nx,Ny);

smart_err_hist = [];

for smartIter = 1:maxOuterIter
    Told = T;

    % recompute material properties with current T
    [Gamma,S] = computeMaterialProperties(XC,YC,T,Nx,Ny);
    [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);

    % face temperatures and gradients for SMART
    [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(T,ge,gw,gn,gs,Nx,Ny);
    [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T,Gamma,bc,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    [dTdx,dTdy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);

    % non-orthogonal correction
    bCorr = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);

    % rebuild upwind system
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        Ee,Ew,En,Es,CE,CW,CN,CS,S,Vc,bCorr,Nx,Ny);
    [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny);
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc);

    % SMART deferred correction using current gradients
    bDC = computeSMARTCorrection(T,Fe,Fw,Fn,Fs,dTdx,dTdy,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    bC = bC + bDC;

    % solve
    T = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T,solver);

    outerErr = max(abs(T(:) - Told(:)));
    smart_err_hist(end+1) = outerErr;
    fprintf('SMART iter %d, max dT = %e, max bDC = %e\n', smartIter, outerErr, max(abs(bDC(:))));

    if outerErr < outerTol
        fprintf('SMART converged in %d iterations\n', smartIter);
        break
    end
end

T_SMART = T;

%% ---- PLOTS ----
figure
contourf(XC, YC, T_upwind, 20, 'LineColor', 'none')
colorbar
colormap(jet)          % change color scheme here\
title('T - upwind scheme (convection-diffusion)')
xlabel('x'); ylabel('y')

figure
contourf(XC, YC, T_SMART, 20, 'LineColor', 'none')
colorbar
colormap(jet)          % change color scheme here
title('T - SMART scheme (convection-diffusion)')
xlabel('x'); ylabel('y')

% difference between the two
figure
contourf(XC, YC, T_SMART - T_upwind, 20, 'LineColor', 'none')
colorbar
colormap(jet)          % change color scheme here
title('T SMART - T upwind (difference)')
xlabel('x'); ylabel('y')

% convergence history
figure
semilogy(upwind_err_hist, 'b-o', 'DisplayName', 'Upwind')
hold on
semilogy(smart_err_hist, 'g-s', 'DisplayName', 'SMART')
xlabel('iteration'); ylabel('max dT')
title('convergence history')
legend; grid on