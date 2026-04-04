clc
clear
close all

%% Grid size
meshSizes = [10 10; 20 20; 40 40; 80 80];

if size(meshSizes,2) == 1
    meshSizes = [meshSizes meshSizes];
elseif size(meshSizes,1) == 1 && numel(meshSizes) > 2
    meshSizes = [meshSizes(:) meshSizes(:)];
end

numMeshes = size(meshSizes,1);
results = struct([]);

%% Boundary conditions
bc.west.type = 'dirichlet';
bc.west.value = 400.0;

bc.east.type = 'dirichlet';
bc.east.value = 300.0;

bc.south.type = 'neumann';
bc.south.value = 0.0;

bc.north.type = 'neumann';
bc.north.value = 0.0;

for meshIdx = 1:numMeshes
    m = meshSizes(meshIdx,1);
    n = meshSizes(meshIdx,2);

    Nx = m-1;
    Ny = n-1;

    fprintf('\n=== Benchmark mesh %d/%d: m = %d, n = %d ===\n', ...
        meshIdx,numMeshes,m,n);

    %% Initial temperature guess
    T = ones(Nx,Ny);

    %% Generate grid
    [X,Y] = Grid_test(m,n);

    %% Geometry
    [XC,YC] = getCentroids(X,Y);

    elemVolume = getCellVolumes(X,Y);

    [Se,Sw,Sn,Ss] = getSurfaceVectors(X,Y,Nx,Ny);

    [Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys] = getFaceCenters(X,Y,Nx,Ny);

    [CE,CW,CN,CS] = getCFVectors(XC,YC,Nx,Ny);

    %% Interpolation factors
    [ge,gw,gn,gs] = getInterpolationFactors( ...
        Nx,Ny,elemVolume,CE,CW,CN,CS,'distance');

    %% Material properties
    %[Gamma,S] = computeMaterialProperties(XC,YC,T,Nx,Ny);
    Gamma = ones(Nx,Ny);
    S = zeros(Nx,Ny);

    %% Gamma interpolation to faces
    [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = ...
        interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);

    %% Temperature interpolation to faces
    [Te,Tw,Tn,Ts] = ...
        interpolateTemperatureToFaces(T,ge,gw,gn,gs,Nx,Ny);

    [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures( ...
        Te,Tw,Tn,Ts,T,Gamma,bc, ...
        XC,YC, ...
        Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys, ...
        Nx,Ny);

    %% Cell gradients
    [dTdx,dTdy] = computeCellGradient( ...
        Te,Tw,Tn,Ts, ...
        Se,Sw,Sn,Ss, ...
        elemVolume,Nx,Ny);

    %% Surface vector decomposition
    [Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = ...
        decomposeSurfaceVectors( ...
        Se,Sw,Sn,Ss, ...
        CE,CW,CN,CS, ...
        Nx,Ny,'overrelaxed');

    %% Non-orthogonal correction
    bCorr = computeNonOrthogonalCorrection( ...
        Gamma_e,Gamma_w,Gamma_n,Gamma_s, ...
        dTdx,dTdy, ...
        Tev,Twv,Tnv,Tsv, ...
        ge,gw,gn,gs, ...
        Nx,Ny);

    %% Diffusion coefficients
    [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients( ...
        Gamma_e,Gamma_w,Gamma_n,Gamma_s, ...
        Ee,Ew,En,Es, ...
        CE,CW,CN,CS, ...
        S,elemVolume,bCorr, ...
        Nx,Ny);

    [aE,aW,aN,aS,aC,bC] = applyAnalyticalBoundaryConditions( ...
        aE,aW,aN,aS,aC,bC, ...
        Gamma, ...
        XC,YC, ...
        Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys, ...
        Se,Sw,Sn,Ss, ...
        Nx,Ny);

    %% Assemble global system
    [A,b] = assembleSystem(aE,aW,aN,aS,aC,bC,Nx,Ny);
    Tvec = A\b;
    Tsol = reshape(Tvec,Nx,Ny);

    residual = A*Tvec - b;
    maxResidual = max(abs(residual));

    fprintf('Max residual = %e\n', maxResidual);

    [L2_error,Max_error,T_exact,errorField] = compareWithAnalytical( ...
        Tsol,X,Y,Nx,Ny,false);

    if numMeshes == 1 || meshIdx == numMeshes
        plotAnalyticalComparison(Tsol,T_exact,errorField,XC,YC,m,n);
    end

    results(meshIdx).m = m;
    results(meshIdx).n = n;
    results(meshIdx).MaxResidual = maxResidual;
    results(meshIdx).L2Error = L2_error;
    results(meshIdx).MaxError = Max_error;
end
if numMeshes > 1
    disp(struct2table(results));
end
