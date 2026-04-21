clc
clear
close all 
m =20; n =20;
Nx = m-1; Ny = n-1; 
u = 10/sqrt(2);
v = 10/sqrt(2);
rho = 1; 
maxOuterIter = 200;
outerTol = 1e-8; 
solver.type = 'TDMA';
solver.maxIter = 5000;
solver.tol = 1e-10; 
  IntFactor = "volume";
Decomp_method = "minimum";

%% BCs
% choose for the west
bc.west.type = 'dirichlet';
bc.west.value = 400.0;
% choose for the east
bc.east.type = 'convective';
bc.east.value = 0.0;
% choose for the south
bc.south.type = 'dirichlet';
bc.south.value = 320.0;
% choose for the north
bc.north.type = 'neumann';
bc.north.value = 0.0;

%Initial conditions
T = 350*ones(Nx,Ny);

%% grid and geometry
[X,Y] = Grid(m,n);
[XC,YC] = getCentroids(X,Y);
Vc = getCellVolumes(X,Y);
[Se,Sw,Sn,Ss]=getSurfaceVectors(X,Y,Nx,Ny);
[Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys]=getFaceCenters(X,Y,Nx,Ny);
[CE,CW,CN,CS]=getCFVectors(XC,YC,Nx,Ny);
[ge,gw,gn,gs]=getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,IntFactor);
[Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv]=decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,Decomp_method);

[Fe,Fw,Fn,Fs] = computeFaceFlux(rho,u,v,Se,Sw,Sn,Ss,Nx,Ny);

upwind_err_hist = [];
for outerIter = 1:maxOuterIter
    Told = T;

    % material properties
    [Gamma,S]=computeMaterialProperties(XC,YC,T,Nx,Ny);

 
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
        [aE,aW,aN,aS,aC,bC]=applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
    XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc,Fe,Fw,Fn,Fs);

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
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny,Fe,Fw,Fn,Fs);
    [dTdx,dTdy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);

    % non-orthogonal correction
    bCorr = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);

    % rebuild upwind system
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        Ee,Ew,En,Es,CE,CW,CN,CS,S,Vc,bCorr,Nx,Ny);
    [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny);
    [aE,aW,aN,aS,aC,bC]=applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
    XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc,Fe,Fw,Fn,Fs);
    bDC = computeSMARTCorrection(T,Fe,Fw,Fn,Fs,dTdx,dTdy,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    bC = bC + bDC;

    % solve
    smartURF = 0.5;
    T_solved = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T,solver);
    T = smartURF*T_solved + (1-smartURF)*T;

    outerErr = max(abs(T(:) - Told(:)));
    smart_err_hist(end+1) = outerErr;
    fprintf('SMART iter %d, max dT = %e, max bDC = %e\n', smartIter, outerErr, max(abs(bDC(:))));

    if outerErr < outerTol
        fprintf('SMART converged in %d iterations\n', smartIter);
        break
    end
end

T_SMART = T;


%% PLEASE NOTE THAT I USED AN LLM TO HELP ME PLOT


%%---- PLOTS ----
figure
contourf(XC, YC, T_upwind, 20, 'LineColor', 'none')
colorbar
colormap(turbo)          % change color scheme here\
title('T - upwind scheme (convection-diffusion)')
xlabel('x'); ylabel('y')
exportgraphics(gcf,'upwind_nonuniform.png','Resolution',300);

figure
contourf(XC, YC, T_SMART, 20, 'LineColor', 'none')
colorbar
colormap(turbo)          % change color scheme here
title('T - SMART scheme (convection-diffusion)')
xlabel('x'); ylabel('y')
exportgraphics(gcf,'smart_nonuniform.png','Resolution',300);

% difference between the two
figure
contourf(XC, YC, T_SMART - T_upwind, 20, 'LineColor', 'none')
colorbar
colormap(turbo)          % change color scheme here
title('T SMART - T upwind (difference)')
xlabel('x'); ylabel('y')
exportgraphics(gcf,'difference_nonuniform.png','Resolution',300);
% convergence history
figure
semilogy(upwind_err_hist, 'b-o', 'DisplayName', 'Upwind')
hold on
semilogy(smart_err_hist, 'g-s', 'DisplayName', 'SMART')
xlabel('iteration'); ylabel('max dT')
title('convergence history')
legend; grid on
exportgraphics(gcf,'convergence_nonuniform.png','Resolution',300);
 