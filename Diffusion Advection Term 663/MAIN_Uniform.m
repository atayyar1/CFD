clc
clear
close all
% grid nodes
m = 31; 
n= 26;
Nx= m-1; 
Ny= n-1;
upwind_err_hist = [];
% velocity, 45 degrees
u= 1/sqrt(2);
v= 1/sqrt(2);
rho= 1;

maxIter = 100;
tol = 1e-3;
solver.type = 'gauss-seidel';
solver.maxIter = 50;
solver.tol = 1e-10;

% BCs
bc.west.type = 'convective'; bc.west.value = 1.0;
bc.north.type = 'convective'; bc.north.value = 1.0;
bc.south.type = 'convective'; bc.south.value = 0.0;
bc.east.type = 'convective';

% initial guess
phi =0.5 * ones(Nx,Ny);

% generate uniform grid
[X,Y] = Grid_test(m,n);

% geometry
[XC,YC] = getCentroids(X,Y);
Vc = getCellVolumes(X,Y);
[Se,Sw,Sn,Ss] = getSurfaceVectors(X,Y,Nx,Ny);
[Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys] = getFaceCenters(X,Y,Nx,Ny);
[CE,CW,CN,CS] = getCFVectors(XC,YC,Nx,Ny);

[ge,gw,gn,gs] = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,"volume");
[Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,"minimum");
[Fe,Fw,Fn,Fs] = computeFaceFlux(rho,u,v,Se,Sw,Sn,Ss,Nx,Ny);
Gamma = zeros(Nx,Ny);
S = zeros(Nx,Ny);
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

    err = max(abs(phi(:) - phi_old(:)));  % ← actual err computed here
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
  bC = bC +1 * bDC;
    
    fprintf('max bDC = %e, max phi change = %e\n', max(abs(bDC(:))), err);
    
    phi = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,phi,solver);
    
    % inside SMART loop, replace:
err = norm(phi(:)-phi_old(:))/norm(phi(:)); 
     smart_err_hist(end+1) = err;

    fprintf('SMART iter %d, err = %e\n', iter, err);
    if err < tol; break; end
end
phi_SMART = phi;
%% contour plots
figure
contourf(XC, YC, phi_upwind, 200, 'LineColor','none')
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

%% PLease note that I used an AI assistant to help me with the plotting functions 



%% ---------------- Uniform case figures ----------------

% Common settings
fs = 12;          % font size
lw = 1.8;         % line width
ms = 6;           % marker size
nLevels = 20;     % contour levels

%% 1) Side-by-side contour plots
figure('Color','w','Position',[100 100 1000 420]);

subplot(1,2,1)
contourf(XC, YC, phi_upwind, nLevels, 'LineColor','none');
colorbar
colormap(parula)
caxis([0 1])
axis equal tight
xlabel('x','FontSize',fs)
ylabel('y','FontSize',fs)
title('Upwind Scheme','FontSize',fs+1)
set(gca,'FontSize',fs,'LineWidth',1.0)
box on

subplot(1,2,2)
contourf(XC, YC, phi_SMART, nLevels, 'LineColor','none');
colorbar
colormap(parula)
caxis([0 1])
axis equal tight
xlabel('x','FontSize',fs)
ylabel('y','FontSize',fs)
colormap(turbo)  
title('SMART Scheme','FontSize',fs+1)
set(gca,'FontSize',fs,'LineWidth',1.0)
box on

sgtitle('Pure Convection Benchmark: Scalar Contours','FontSize',fs+2)

exportgraphics(gcf,'uniform_contours.png','Resolution',300);


%% 2) Centerline comparison
figure('Color','w','Position',[150 150 700 500]);

plot(YC(ic,:), phi_upwind(ic,:), '-o', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerIndices', 1:2:length(YC(ic,:)), ...
    'DisplayName', 'Upwind');
hold on

plot(YC(ic,:), phi_SMART(ic,:), '-s', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerIndices', 1:2:length(YC(ic,:)),'color',  [0.85 0.33 0.10],...
    'DisplayName', 'SMART');

plot(y_analytical, phi_analytical, '--', ...
    'LineWidth', 2.2, 'Color', [0.2 0.2 0.2],'DisplayName', 'Analytical');

xlabel('y','FontSize',fs)
ylabel('\phi','FontSize',fs)
title('Vertical Centerline Profile','FontSize',fs+1)
legend('Location','best','FontSize',fs-1)
box on
xlim([0 1])
ylim([-0.05 1.05])
set(gca,'FontSize',fs,'LineWidth',1.0)

exportgraphics(gcf,'centerline_uniform.png','Resolution',300);


%% 3) Convergence history
figure('Color','w','Position',[200 200 700 500]);

semilogy(upwind_err_hist, '-o', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerIndices', 1:max(1,round(length(upwind_err_hist)/12)):length(upwind_err_hist), ...
    'DisplayName', 'Upwind');
hold on

semilogy(smart_err_hist, '-s', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerIndices', 1:max(1,round(length(smart_err_hist)/12)):length(smart_err_hist), ...
    'DisplayName', 'SMART');

xlabel('Iteration','FontSize',fs)
ylabel('Error','FontSize',fs)
title('Convergence History','FontSize',fs+1)
legend('Location','best','FontSize',fs-1)
box on
set(gca,'FontSize',fs,'LineWidth',1.0)

exportgraphics(gcf,'convergence_uniform.png','Resolution',300);