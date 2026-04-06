clc
clear
close all

%% Grid size
%nx and ny
meshSizes = [20,20];

if size(meshSizes,2) == 1
    meshSizes = [meshSizes meshSizes];
elseif size(meshSizes,1) == 1 && numel(meshSizes) > 2
    meshSizes = [meshSizes(:) meshSizes(:)];
end

numMeshes = size(meshSizes,1);
results = struct([]);
plotData = struct([]);

%% Solver settings
maxOuterIter = 100;
outerTol = 1e-6;

%% Interpolation Factor Method, choose between distance and volume:
IntFactor = "volume";
%% Choose decomposition Method, choose between overrelaxed, minimum, or normal
Decomp_method = "minimum";
%% PLEASE choose between TDMA and Gauss_seidel
solver.type = 'gauss-seidel'; 
solver.maxIter = 5000;
solver.tol = 1e-10;

%% Boundary conditions
% choose for the west
bc.west.type = 'dirichlet';
bc.west.value = 400.0;


% choose for the east
bc.east.type = 'dirichlet';
bc.east.value = 400.0;

% choose for the south
bc.south.type = 'dirichlet';
bc.south.value = 400;

% choose for the north
bc.north.type = 'dirichlet';
bc.north.value = 400;

for meshIdx = 1:numMeshes
    m = meshSizes(meshIdx,1);
    n = meshSizes(meshIdx,2);

    Nx = m-1;
    Ny = n-1;

    fprintf('\n=== Mesh %d/%d: m = %d, n = %d ===\n', ...
        meshIdx,numMeshes,m,n);

    %% Initial temperature guess
    T = 350*ones(Nx,Ny);

    %% Generate grid
    [X,Y] = Grid(m,n);

    %% Geometry
    [XC,YC] = getCentroids(X,Y);

    elemVolume = getCellVolumes(X,Y);

    [Se,Sw,Sn,Ss] = getSurfaceVectors(X,Y,Nx,Ny);

    [Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys] = getFaceCenters(X,Y,Nx,Ny);

    [CE,CW,CN,CS] = getCFVectors(XC,YC,Nx,Ny);

    %% Interpolation factors
    [ge,gw,gn,gs] = getInterpolationFactors( ...
        Nx,Ny,elemVolume,CE,CW,CN,CS,IntFactor);

    %% Surface vector decomposition
    [Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = ...
        decomposeSurfaceVectors( ...
        Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,Decomp_method);

    %% Nonlinear iterations
    for outerIter = 1:maxOuterIter

        Told = T;

        %% Material properties
        [Gamma,S] = computeMaterialProperties(XC,YC,T,Nx,Ny);

        %% Gamma interpolation to faces
        [Gamma_e,Gamma_w,Gamma_n,Gamma_s]= ...
            interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny);

        %% Temperature interpolation to faces
        [Te,Tw,Tn,Ts] = ...
            interpolateTemperatureToFaces(T,ge,gw,gn,gs,Nx,Ny);

        [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures( ...
            Te,Tw,Tn,Ts,T,Gamma,bc, ...
            XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny);

        %% Cell gradients
        [dTdx,dTdy] = computeCellGradient( ...
            Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,elemVolume,Nx,Ny);

        %% Non-orthogonal correction
        % I have other ones because it caused me so many problems, so I started
        % testing, I am keeping them for referencing
        bCorr = computeNonOrthogonalCorrection( ...
            Gamma_e,Gamma_w,Gamma_n,Gamma_s,dTdx,dTdy, ...
            Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny);
        %bCorr =zeros(Nx,Ny); % Disable non-orthogonal correction for testing
        %bCorr = bCorr*0.5; % Relax the non-orthogonal correction to improve convergence

        %% Diffusion coefficients
        [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients( ...
            Gamma_e,Gamma_w,Gamma_n,Gamma_s,Ee,Ew,En,Es, ...
            CE,CW,CN,CS,S,elemVolume,bCorr,Nx,Ny);

        [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions( ...
            aE,aW,aN,aS,aC,bC,  Gamma, XC,YC, ...
            Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc);

        %% Linear solve
        T = solveLinearSystem(aE,aW,aN,aS,aC,bC,Nx,Ny,T,solver);

        outerErr = max(abs(T(:) - Told(:)));
        fprintf('Nonlinear iteration %d, max dT = %e\n', outerIter, outerErr);

        if outerErr < outerTol
            fprintf('Nonlinear iterations converged in %d iterations\n', outerIter);
            break
        end

    end

    Tsol = T;

    %% Assemble final system for residual check
    [A,b] = assembleSystem(aE,aW,aN,aS,aC,bC,Nx,Ny);
    Tvec = Tsol(:);

    % I had problems with the residuals being very large, so I added this to check if there is a problem with the assembly or the solve, 
    % and it seems that the residuals are reasonable, which is good news
   
    residual = A*Tvec - b;
    maxResidual = max(abs(residual));
    Tmin = min(Tsol(:));
    Tmax = max(Tsol(:));


    %%Please note that I used codex for the plots and the results display, so I have not added comments to them, if you want me to add comments please let me know.
    fprintf('Max residual = %e\n', maxResidual);
    fprintf('Temperature range = [%f, %f]\n', Tmin, Tmax);

    results(meshIdx).m = m;
    results(meshIdx).n = n;
    results(meshIdx).OuterIterations = outerIter;
    results(meshIdx).MaxResidual = maxResidual;
    results(meshIdx).Tmin = Tmin;
    results(meshIdx).Tmax = Tmax;

    plotData(meshIdx).m = m;
    plotData(meshIdx).n = n;
    plotData(meshIdx).XC = XC;
    plotData(meshIdx).YC = YC;
    plotData(meshIdx).Tsol = Tsol;
end

if numMeshes > 1
    % Display results in a table, converting the struct array to a table for better formatting
    disp(struct2table(results));
end

plotTemperatureMeshes(plotData,IntFactor,Decomp_method,solver.type);
% plot the vectors
