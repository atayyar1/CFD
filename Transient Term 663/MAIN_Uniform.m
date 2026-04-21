clc
clear
close all

%% parameters
a_dimension = 0.1;
b_dimension = 0.2
m = 21; n = 41
Nx = m-1; Ny = n-1;

rho    = 7800;
cp     = 500;
dt     = 1;
t_end  = 60;
Nsteps = t_end/dt;

save_times = [10 30 60];

solver.type    = 'TDMA'; % for the user to choose, TDMA or GAUSS SEIDEL
solver.maxIter = 5000;
solver.tol     = 1e-10;

IntFactor    = "volume";
Decomp_method = "minimum";

%% BCs
bc.west.type  = 'neumann';   bc.west.value  = 0.0;
bc.east.type  = 'neumann';   bc.east.value  = 0.0;
bc.south.type = 'dirichlet'; bc.south.value = 0.0;
bc.north.type = 'dirichlet'; bc.north.value = 0.0;

%% grid and geometry
[X,Y]    = GridRect(m, n, a_dimension, b_dimension);
[XC,YC]  = getCentroids(X,Y);
Vc     = getCellVolumes(X,Y);
[Se,Sw,Sn,Ss]                     = getSurfaceVectors(X,Y,Nx,Ny);
[Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys]        = getFaceCenters(X,Y,Nx,Ny);
[CE,CW,CN_vec,CS]                 = getCFVectors(XC,YC,Nx,Ny);
[ge,gw,gn,gs]                     = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN_vec,CS,IntFactor);
[Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv]    = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN_vec,CS,Nx,Ny,Decomp_method);

% constant gamma for this problem, later, this will changed to the one in
% assignment 1, on the non-orthogonal grid, this will be done for the test
% case on the uniform grid
Gamma = 15*ones(Nx,Ny);

% storage
T_be_save   = zeros(Nx,Ny,length(save_times));
T_cn_save   = zeros(Nx,Ny,length(save_times));
T_anal_save = zeros(Nx,Ny,length(save_times));
T_be = 1000 * ones(Nx,Ny);

for step = 1:Nsteps
    t    = step * dt;
    Told = T_be;

    % face gammas (constant here but keeping it general)
    [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);

    % face temps and cell gradients for non-ortho correction
    [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(T_be,ge,gw,gn,gs,Nx,Ny);
    [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T_be,Gamma,bc,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    [dTdx,dTdy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);
    bCorr      = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);

    % diffusion system - S=0 for this problem, also, same as Gamma
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        Ee,Ew,En,Es,CE,CW,CN_vec,CS,zeros(Nx,Ny),Vc,bCorr,Nx,Ny);

    % BCs - no convection so dont pass face fluxes
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc);

    %transient term
    [aC,bC] = assembleTransientTerm(Told,Vc,rho,cp,dt,aC,bC,Nx,Ny);

    % solve
    T_be = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T_be,solver);

    % save at required times
    for kk = 1:length(save_times)
        if t == save_times(kk)
            T_be_save(:,:,kk)   = T_be;
            T_anal_save(:,:,kk) = computeAnalyticalSolution(XC,YC,t,a_dimension,b_dimension,15,rho,cp);
        end
    end
end

%% Now, I will add Crank Nicolson
fprintf('Crank Nicholson\n')
  T_cn      = 1000 * ones(Nx,Ny);
T_be_prev = 1000 * ones(Nx,Ny);   % need two time levels back

for step = 1:Nsteps
    t    = step * dt;
   Told = T_be_prev;


    % face gammas
    [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);

    % face temps and gradients for non-ortho correction
    [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(T_be_prev,ge,gw,gn,gs,Nx,Ny);
[Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T_be_prev,Gamma,bc,...
    XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);
    [dTdx,dTdy]   = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny);
    bCorr         = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);

    % build diffusion system
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,...
        Ee,Ew,En,Es,CE,CW,CN_vec,CS,zeros(Nx,Ny),Vc,bCorr,Nx,Ny);
    [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,...
        XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc);

    % add transient - same as BE (step 1)
    [aC,bC] = assembleTransientTerm(Told,Vc,rho,cp,dt,aC,bC,Nx,Ny);

    % step 1: solve BE
T_be = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T_be_prev,solver);
    % step 2: CN extrapolation - no second solve needed
    % T_cn^{t+dt} = 2*T_be^t - T_be^{t-dt}
    T_cn = 2*T_be - T_be_prev;

    % update previous BE level for next step
    T_be_prev = T_be;
1
    % save at required times
    for kk = 1:length(save_times)
        if t == save_times(kk)
            T_cn_save(:,:,kk) = T_cn;
            fprintf('  CN: saved t = %d s\n', t)
        end
    end
end


%% PLEASE NOTE: I have used an LLM To help me with the ploting functions! 
%% ---- PLOTS ----
plotTransientResults(XC, YC, T_be_save, T_cn_save, T_anal_save, save_times);