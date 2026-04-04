clc
clear
close all

% grid nodes
m = 31; n = 26;
Nx = m-1; Ny = n-1;
upwind_err_hist = [];

% velocity field, 45 degrees
u = 1/sqrt(2);
v = 1/sqrt(2);
rho = 1;

% solver settings
maxIter = 200;
tol = 1e-8;

solver.type = 'gauss-seidel';
solver.maxIter = 5000;
solver.tol = 1e-10;

% BCs
bc.west.type = 'convective'; bc.west.value = 1.0;
bc.north.type = 'convective'; bc.north.value = 1.0;
bc.south.type = 'convective'; bc.south.value = 0.0;
bc.east.type = 'convective';

% initial guess
phi = zeros(Nx,Ny);

% generate uniform grid
[X,Y] = Grid_test(m,n);

% geometry
[XC,YC] = getCentroids(X,Y);
Vc = getCellVolumes(X,Y);
[Se,Sw,Sn,Ss] = getSurfaceVectors(X,Y,Nx,Ny);
[Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys] = getFaceCenters(X,Y,Nx,Ny);
[CE,CW,CN,CS] = getCFVectors(XC,YC,Nx,Ny);

% interpolation factors
[ge,gw,gn,gs] = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,"volume");

% surface vector decomposition (needed for function calls even if correction~=0)
[Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,"minimum");
% face fluxes, constant so computed once outside the loop
[Fe,Fw,Fn,Fs] = computeFaceFlux(rho,u,v,Se,Sw,Sn,Ss,Nx,Ny);
% pure convection, Gamma = 0
Gamma = zeros(Nx,Ny);
S = zeros(Nx,Ny);
%% main loop
for iter = 1:maxIter
    phi_old = phi;
    % diffusion coefficients (all zeros since Gamma=0, just initializes arrays)
    bCorr = zeros(Nx,Ny);
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma,Gamma,Gamma,Gamma, ...
        Ee,Ew,En,Es,CE,CW,CN,CS,S,Vc,bCorr,Nx,Ny);

    % add convection
    [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny);
 
    % boundary conditions
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc,Fe,Fw,Fn,Fs);
 
    % solve
    phi = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,phi,solver);

    err = max(abs(phi(:) - phi_old(:)));
    fprintf('iter %d, err = %e\n', iter, err);
    upwind_err_hist(end+1) = err;
    if err < tol
        fprintf('converged at iteration %d\n', iter);
        break
    end

end

phi_upwind = phi;
smart_err_hist = [];

% SMART loop
for iter = 1:maxIter
    phi_old = phi;
    
    % recompute gradients using current phi
    [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(phi,ge,gw,gn,gs,Nx,Ny);
    [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,phi,Gamma,bc,...
    XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    [dphidx,dphidy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);
    
    % same upwind matrix
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma,Gamma,Gamma,Gamma,...
        Ee,Ew,En,Es,CE,CW,CN,CS,S,Vc,zeros(Nx,Ny),Nx,Ny);
    [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny);
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,XC,YC,...
        Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc,Fe,Fw,Fn,Fs);
    
    % add SMART correction to RHS
    bDC = computeSMARTCorrection(phi,Fe,Fw,Fn,Fs,dphidx,dphidy,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    fprintf('SMART iter %d, max correction = %e, err = %e\n', iter, max(abs(bDC(:))), err);
    bC = bC +  0.95*bDC;
    
    phi = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,phi,solver);
    
    % inside SMART loop, replace:
err = max(abs(phi(:) - phi_old(:))); 
     smart_err_hist(end+1) = err;

    fprintf('SMART iter %d, err = %e\n', iter, err);
    if err < tol; break; end
end

phi_SMART = phi;
%% contour plots
figure
contourf(XC, YC, phi_upwind, 20, 'LineColor','none')
colorbar
title('phi - upwind scheme')
xlabel('x'); ylabel('y')

figure
contourf(XC, YC, phi_SMART, 20, 'LineColor','none')
colorbar
title('phi - SMART scheme')
xlabel('x'); ylabel('y')

%% centerline comparison
[~, ic] = min(abs(XC(:,1) - 0.5));

y_analytical = [0, 0.5, 0.5, 1];
phi_analytical = [0, 0, 1, 1];

figure
plot(YC(ic,:), phi_upwind(ic,:), 'b-o', 'DisplayName', 'Upwind')
hold on
plot(YC(ic,:), phi_SMART(ic,:), 'g-s', 'DisplayName', 'SMART')
plot(y_analytical, phi_analytical, 'r--', 'DisplayName', 'Analytical')
xlabel('y'); ylabel('\phi')
title('centerline profile - upwind vs SMART vs analytical')
legend; grid on

%% Plot convergence
figure
semilogy(upwind_err_hist, 'b-o', 'DisplayName', 'Upwind')
hold on
% your current plot should already use semilogy
% if it's showing linear scale, check your plot command
semilogy(smart_err_hist, 'g-s', 'DisplayName', 'SMART')
xlabel('iteration'); ylabel('max error')
title('convergence history')
legend; grid on
