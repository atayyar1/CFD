clc
clear
close all

%% parameters
m = 21; n = 21;
Nx = m-1; Ny = n-1;

u = 10/sqrt(2);
v = 10/sqrt(2);
rho    = 1;
cp     = 500;
dt     = 1;
t_end  = 200;
Nsteps = t_end/dt;

save_times = [10 30 200];

solver.type    = 'TDMA';
solver.maxIter = 5000;
solver.tol     = 1e-10;

IntFactor     = "volume";
Decomp_method = "minimum";

%% BCs — A1 domain
bc.west.type  = 'dirichlet'; bc.west.value  = 400.0;
bc.east.type  = 'convective'; bc.east.value  = 0.0;
bc.south.type = 'dirichlet'; bc.south.value = 320.0;
bc.north.type = 'neumann';   bc.north.value = 0.0;

%% grid and geometry
[X,Y]    = Grid(m,n);
[XC,YC]  = getCentroids(X,Y);
Vc       = getCellVolumes(X,Y);
[Se,Sw,Sn,Ss]                  = getSurfaceVectors(X,Y,Nx,Ny);
[Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys]     = getFaceCenters(X,Y,Nx,Ny);
[CE,CW,CN,CS]                  = getCFVectors(XC,YC,Nx,Ny);
[ge,gw,gn,gs]                  = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,IntFactor);
[Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,Decomp_method);

%% face fluxes — constant velocity so computed once outside loop
[Fe,Fw,Fn,Fs] = computeFaceFlux(rho,u,v,Se,Sw,Sn,Ss,Nx,Ny);

%% storage
T_be_save = zeros(Nx,Ny,length(save_times));
T_cn_save = zeros(Nx,Ny,length(save_times));
time_labels = {'t = 10s', 't = 30s', 't = 60s'};

%% ---- BACKWARD EULER ----
fprintf('Backward Euler\n')
T_be = 1000 * ones(Nx,Ny);

for step = 1:Nsteps
    t    = step * dt;
    Told = T_be;

    % variable material properties
    [Gamma, S] = computeMaterialProperties(XC, YC, T_be, Nx, Ny);

    % face gammas
    [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);

    % face temperatures for non-ortho correction and SMART
    [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(T_be,ge,gw,gn,gs,Nx,Ny);
    [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T_be,Gamma,bc,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny,Fe,Fw,Fn,Fs);

    % cell gradients
    [dTdx,dTdy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);

    % non-ortho correction
    bCorr = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);

    % diffusion system
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        Ee,Ew,En,Es,CE,CW,CN,CS,S,Vc,bCorr,Nx,Ny);

    % upwind convection
    [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny);

    % BCs
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc,Fe,Fw,Fn,Fs);

    % SMART deferred correction
    bDC = computeSMARTCorrection(T_be,Fe,Fw,Fn,Fs,dTdx,dTdy,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    bC = bC + bDC;

    % transient term
    [aC,bC] = assembleTransientTerm(Told,Vc,rho,cp,dt,aC,bC,Nx,Ny);

    % solve
    T_be = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T_be,solver);

    % save
    for kk = 1:length(save_times)
        if t == save_times(kk)
            T_be_save(:,:,kk) = T_be;
            fprintf('  BE: saved t = %d s\n', t)
        end
    end
end

%% ---- CRANK-NICOLSON ----
fprintf('Crank Nicolson\n')
T_cn      = 1000 * ones(Nx,Ny);
T_be_prev = 1000 * ones(Nx,Ny);

for step = 1:Nsteps
    t    = step * dt;
    Told = T_be_prev;   % T_be_prev, NOT T_cn

    % material properties based on previous BE level
    [Gamma, S] = computeMaterialProperties(XC, YC, T_be_prev, Nx, Ny);

    % face gammas
    [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);

    % face temperatures — use T_be_prev as the field for this BE step
    [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(T_be_prev,ge,gw,gn,gs,Nx,Ny);
    [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T_be_prev,Gamma,bc,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny,Fe,Fw,Fn,Fs);

    % cell gradients
    [dTdx,dTdy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);

    % non-ortho correction
    bCorr = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);

    % diffusion system
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        Ee,Ew,En,Es,CE,CW,CN,CS,S,Vc,bCorr,Nx,Ny);

    % upwind convection
    [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny);

    % BCs
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc,Fe,Fw,Fn,Fs);

    % SMART deferred correction
    bDC = computeSMARTCorrection(T_be_prev,Fe,Fw,Fn,Fs,dTdx,dTdy,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    bC = bC + bDC;

    % transient term
    [aC,bC] = assembleTransientTerm(Told,Vc,rho,cp,dt,aC,bC,Nx,Ny);

    % step 1: BE solve
    T_be = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T_be_prev,solver);

    % step 2: CN extrapolation — no second solve
    T_cn = 2*T_be - T_be_prev;

    % update
    T_be_prev = T_be;

    % save
    for kk = 1:length(save_times)
        if t == save_times(kk)
            T_cn_save(:,:,kk) = T_cn;
            fprintf('  CN: saved t = %d s\n', t)
        end
    end
end

%% ---- PLOTS ----
%% PLEASE NOTE: I have used an LLM to help me with the plotting functions!
plotNonUniformResults(XC, YC, T_be_save, T_cn_save, save_times, 'figures');