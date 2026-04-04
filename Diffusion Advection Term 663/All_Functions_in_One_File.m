%%All functions:
function [XC, YC] = getCentroids(X, Y)

[m,n] = size(X);

XC = zeros(m-1,n-1);
YC = zeros(m-1,n-1);

for i = 1:m-1
    for j = 1:n-1

        XC(i,j) = 0.25*( X(i,j)   + X(i+1,j) + X(i,j+1) + X(i+1,j+1) );

        YC(i,j) = 0.25*( Y(i,j)   + Y(i+1,j) + Y(i,j+1) + Y(i+1,j+1) );

    end
end

end
function [Gamma, S] = computeMaterialProperties(XC,YC,T,Nx,Ny)
Gamma = zeros(Nx,Ny);
S     = zeros(Nx,Ny);

%Loop over control volumes
for i = 1:Nx
    for j = 1:Ny

        x = XC(i,j);
        y = YC(i,j);

        % Source term
        S(i,j) = T(i,j)*( 2*x - 0.2*y )/400;
        % Diffusion coefficients
        Gamma(i,j) = T(i,j)*( x^2 + exp(0.1*y) )/400;


    end
end

end
function [Te,Tw,Tn,Ts] = interpolateTemperatureToFaces(T,ge,gw,gn,gs,Nx,Ny)
Te = zeros(Nx,Ny);
Tw = zeros(Nx,Ny);
Tn = zeros(Nx,Ny);
Ts = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        %% East
        if i < Nx
            Te(i,j) = ge(i,j)*T(i,j) + (1-ge(i,j))*T(i+1,j);
        end
        %% West
        if i > 1
            Tw(i,j) = gw(i,j)*T(i,j) + (1-gw(i,j))*T(i-1,j);
        end
        %% North
        if j < Ny
            Tn(i,j) = gn(i,j)*T(i,j) + (1-gn(i,j))*T(i,j+1);
        end
        %% South
        if j > 1
            Ts(i,j) = gs(i,j)*T(i,j) + (1-gs(i,j))*T(i,j-1);
        end
    end
end
end
function elemVolume = getCellVolumes(X,Y)

[m,n] = size(X);

Nx = m-1;
Ny = n-1;
elemVolume = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny

        % Nodes of the cell
        x1 = X(i,j);
        y1 = Y(i,j);

        x2 = X(i+1,j);
        y2 = Y(i+1,j);

        x3 = X(i+1,j+1);
        y3 = Y(i+1,j+1);

        x4 = X(i,j+1);
        y4 = Y(i,j+1);

        % Triangle 1
        A1 = 0.5 * abs((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
        % Triangle 2
        A2 = 0.5 * abs((x3-x1)*(y4-y1) - (y3-y1)*(x4-x1));

        % Cell area
        elemVolume(i,j) = A1 + A2;

    end
end

end
function [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny)

Gamma_e = zeros(Nx,Ny);
Gamma_w = zeros(Nx,Ny);
Gamma_n = zeros(Nx,Ny);
Gamma_s = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        if i < Nx
            Gamma_e(i,j)=1/((1-ge(i,j))/Gamma(i+1,j)+ge(i,j)/Gamma(i,j));
        end

        if i > 1
            Gamma_w(i,j)=1/((1-gw(i,j))/Gamma(i-1,j)+gw(i,j)/Gamma(i,j));
        end

        if j < Ny
            Gamma_n(i,j)=1/((1-gn(i,j))/Gamma(i,j+1)+gn(i,j)/Gamma(i,j));
        end

        if j > 1
            Gamma_s(i,j)=1/( (1-gs(i,j))/Gamma(i,j-1)+gs(i,j)/Gamma(i,j));
        end
    end
end
end
function [CE,CW,CN,CS] = getCFVectors(XC,YC,Nx,Ny)

CE = zeros(Nx,Ny,2);
CW = zeros(Nx,Ny,2);
CN = zeros(Nx,Ny,2);
CS = zeros(Nx,Ny,2);

for i = 1:Nx
    for j = 1:Ny

        %% East
        if i < Nx
            CE(i,j,1) = XC(i+1,j) - XC(i,j);
            CE(i,j,2) = YC(i+1,j) - YC(i,j);
        end

        %% WEst
        if i > 1
            CW(i,j,1) = XC(i-1,j) - XC(i,j);
            CW(i,j,2) = YC(i-1,j) - YC(i,j);
        end

        %% North
        if j < Ny
            CN(i,j,1) = XC(i,j+1) - XC(i,j);
            CN(i,j,2) = YC(i,j+1) - YC(i,j);
        end

        %% SOUTH
        if j > 1
            CS(i,j,1) = XC(i,j-1) - XC(i,j);
            CS(i,j,2) = YC(i,j-1) - YC(i,j);
        end

    end
end

end
function [Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys] = getFaceCenters(X,Y,Nx,Ny)

Xe = zeros(Nx,Ny);
Ye = zeros(Nx,Ny);

Xw = zeros(Nx,Ny);
Yw = zeros(Nx,Ny);

Xn = zeros(Nx,Ny);
Yn = zeros(Nx,Ny);

Xs = zeros(Nx,Ny);
Ys = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny

        %% EAST FACE
        Xe(i,j) = 0.5*( X(i+1,j) + X(i+1,j+1) );
        Ye(i,j) = 0.5*( Y(i+1,j) + Y(i+1,j+1) );

        %% WEST FACE
        Xw(i,j) = 0.5*( X(i,j) + X(i,j+1) );
        Yw(i,j) = 0.5*( Y(i,j) + Y(i,j+1) );

        %% NORTH FACE
        Xn(i,j) = 0.5*( X(i,j+1) + X(i+1,j+1) );
        Yn(i,j) = 0.5*( Y(i,j+1) + Y(i+1,j+1) );

        %% SOUTH FACE
        Xs(i,j) = 0.5*( X(i,j) + X(i+1,j) );
        Ys(i,j) = 0.5*( Y(i,j) + Y(i+1,j) );

    end
end

end
function [Se,Sw,Sn,Ss] = getSurfaceVectors(X,Y,Nx,Ny)

Se = zeros(Nx,Ny,2);
Sw = zeros(Nx,Ny,2);
Sn = zeros(Nx,Ny,2);
Ss = zeros(Nx,Ny,2);

for i = 1:Nx
    for j = 1:Ny

        %% EAST FACE
        dx = X(i+1,j+1) - X(i+1,j);
        dy = Y(i+1,j+1) - Y(i+1,j);

        Se(i,j,1) = dy;
        Se(i,j,2) = -dx;

        %% WEST FACE
        dx = X(i,j+1) - X(i,j);
        dy = Y(i,j+1) - Y(i,j);

        Sw(i,j,1) = -dy;
        Sw(i,j,2) =  dx;

        %% NORTH FACE
        dx = X(i+1,j+1) - X(i,j+1);
        dy = Y(i+1,j+1) - Y(i,j+1);

        Sn(i,j,1) = -dy;
        Sn(i,j,2) =  dx;
        %% SOUTH FACE
        dx = X(i+1,j) - X(i,j);
        dy = Y(i+1,j) - Y(i,j);

        Ss(i,j,1) = dy;
        Ss(i,j,2) = -dx;

    end
end

end
function [X,Y] = Grid(m,n)

xi  = linspace(0,1,m);
eta = linspace(0,1,n);

X = zeros(m,n);
Y = zeros(m,n);

for i = 1:m
    Xi = xi(i);

    for j = 1:n
        Eta = eta(j);

        XY = (1-Eta)*Xb(Xi) + Eta*Xt(Xi) + (1-Xi)*Xl(Eta) + Xi*Xr(Eta) - (Xi*Eta*Xt(1) + Xi*(1-Eta)*Xb(1) + ...
            Eta*(1-Xi)*Xt(0) + (1-Xi)*(1-Eta)*Xb(0));

        X(i,j) = XY(1);
        Y(i,j) = XY(2);

    end
end

end

function p = Xb(s)

x = 0 + s*(2-0);
y = 0 + s*(0-0);

p = [x; y];

end


% Top:B-C
function p = Xt(s)

x = -1 + s*(3+1);
y = 3  + s*(5-3);

p = [x; y];

end


% Left: AB
function p = Xl(s)

x = 0  + s*(-1-0);
y = 0  + s*(3-0);

p = [x; y];

end


% Right: DC
function p = Xr(s)

x = 2 + s*(3-2);
y = 0 + s*(5-0);

p = [x; y];

end
function plotCFVectors(X,Y,XC,YC,CE,CW,CN,CS,Nx,Ny)

figure
hold on
axis equal
box on

%% Plot grid
[m,n] = size(X);

for i = 1:m
    plot(X(i,:),Y(i,:),'k')
end

for j = 1:n
    plot(X(:,j),Y(:,j),'k')
end

%% Plot centroids
plot(XC(:),YC(:),'bo','MarkerSize',6,'LineWidth',1.5)

%% Plot CF vectors
scale = 0.4;

for i = 1:Nx
    for j = 1:Ny

        xc = XC(i,j);
        yc = YC(i,j);

        if i < Nx
            quiver(xc,yc,scale*CE(i,j,1),scale*CE(i,j,2),'r')
        end

        if i > 1
            quiver(xc,yc,scale*CW(i,j,1),scale*CW(i,j,2),'g')
        end

        if j < Ny
            quiver(xc,yc,scale*CN(i,j,1),scale*CN(i,j,2),'b')
        end

        if j > 1
            quiver(xc,yc,scale*CS(i,j,1),scale*CS(i,j,2),'m')
        end

    end
end

title('Center-to-Center Vectors (CF)')
xlabel('x')
ylabel('y')

hold off

end
function [Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,method)

Ee  = zeros(Nx,Ny,2);
Ew  = zeros(Nx,Ny,2);
En  = zeros(Nx,Ny,2);
Es  = zeros(Nx,Ny,2);

Tev = zeros(Nx,Ny,2);
Twv = zeros(Nx,Ny,2);
Tnv = zeros(Nx,Ny,2);
Tsv = zeros(Nx,Ny,2);

for i = 1:Nx
    for j = 1:Ny

        %% East
        if i < Nx
            sf = [Se(i,j,1) Se(i,j,2)];
            cf = [CE(i,j,1) CE(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            Ee(i,j,:)  = Ef;
            Tev(i,j,:) = Tf;
        end

        %% West
        if i > 1
            sf = [Sw(i,j,1) Sw(i,j,2)];
            cf = [CW(i,j,1) CW(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            Ew(i,j,:)  = Ef;
            Twv(i,j,:) = Tf;
        end

        %% North
        if j < Ny
            sf = [Sn(i,j,1) Sn(i,j,2)];
            cf = [CN(i,j,1) CN(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            En(i,j,:)  = Ef;
            Tnv(i,j,:) = Tf;
        end

        %% South
        if j > 1
            sf = [Ss(i,j,1) Ss(i,j,2)];
            cf = [CS(i,j,1) CS(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            Es(i,j,:)  = Ef;
            Tsv(i,j,:) = Tf;
        end

    end
end

end

function [Ef,Tf] = oneFaceDecomposition(sf,cf,method)

ef = cf / norm(cf);              % unit vector along CF
Sfmag = norm(sf);
n = sf / Sfmag;                  % unit normal direction of the face
cosTheta = dot(n,ef);

if strcmp(method,'minimum')
    Ef = dot(sf,ef) * ef;

elseif strcmp(method,'normal')
    Ef = Sfmag * ef;

elseif strcmp(method,'overrelaxed')
    Ef = (Sfmag / cosTheta) * ef;

else
    error('Unknown method. Use minimum, normal, or overrelaxed.');
end

Tf = sf - Ef;

end

function [ge,gw,gn,gs] = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,method)


% Initialize arrays
ge = zeros(Nx,Ny);
gw = zeros(Nx,Ny);
gn = zeros(Nx,Ny);
gs = zeros(Nx,Ny);

% Loop over elements
if strcmp(method,'volume')

    % volume-based interpolation
    for i = 1:Nx
        for j = 1:Ny

            if i < Nx
                ge(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i+1,j));
            end

            if i > 1
                gw(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i-1,j));
            end

            if j < Ny
                gn(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i,j+1));
            end

            if j > 1
                gs(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i,j-1));
            end

        end
    end


elseif strcmp(method,'distance')

    % distance-based interpolation
    for i = 1:Nx
        for j = 1:Ny

            if i < Nx
                dP = norm(CE(i,j,:));
                dE = norm(CW(i+1,j,:));
                ge(i,j) = dE/(dP + dE);
            end

            if i > 1
                dP = norm(CW(i,j,:));
                dW = norm(CE(i-1,j,:));
                gw(i,j) = dW/(dP + dW);
            end

            if j < Ny
                dP = norm(CN(i,j,:));
                dN = norm(CS(i,j+1,:));
                gn(i,j) = dN/(dP + dN);
            end

            if j > 1
                dP = norm(CS(i,j,:));
                dS = norm(CN(i,j-1,:));
                gs(i,j) = dS/(dP + dS);
            end

        end
    end

end
end
    function bCorr = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s, ...
            dTdx,dTdy,Tev,Twv,Tnv,Tsv,Nx,Ny)

        bCorr = zeros(Nx,Ny);

        for i = 1:Nx
            for j = 1:Ny

                corr = 0;

                %% East fac
                if i < Nx
                    gradTf = dTdx(i,j)*Tev(i,j,1) + dTdy(i,j)*Tev(i,j,2);
                    corr = corr + Gamma_e(i,j)*gradTf;
                end

                %% West face
                if i > 1
                    gradTf = dTdx(i,j)*Twv(i,j,1) + dTdy(i,j)*Twv(i,j,2);
                    corr = corr + Gamma_w(i,j)*gradTf;
                end

                %% North face
                if j < Ny
                    gradTf = dTdx(i,j)*Tnv(i,j,1) + dTdy(i,j)*Tnv(i,j,2);
                    corr = corr + Gamma_n(i,j)*gradTf;
                end

                %% South face
                if j > 1
                    gradTf = dTdx(i,j)*Tsv(i,j,1) + dTdy(i,j)*Tsv(i,j,2);
                    corr = corr + Gamma_s(i,j)*gradTf;
                end

                bCorr(i,j) = corr;

            end
        end

    end
    function [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,Ee,Ew,En,Es,CE,CW,CN,CS, ...
            S,Vc,bCorr,Nx,Ny)

        aE = zeros(Nx,Ny);
        aW = zeros(Nx,Ny);
        aN = zeros(Nx,Ny);
        aS = zeros(Nx,Ny);
        aC = zeros(Nx,Ny);
        bC = zeros(Nx,Ny);

        for i = 1:Nx
            for j = 1:Ny
                % East
                if i < Nx
                    aE(i,j) = -Gamma_e(i,j) * norm(squeeze(Ee(i,j,:))) / norm(squeeze(CE(i,j,:)));
                end
                %West
                if i > 1
                    aW(i,j) = -Gamma_w(i,j) * norm(squeeze(Ew(i,j,:))) / norm(squeeze(CW(i,j,:)));
                end

                %North
                if j < Ny
                    aN(i,j) = -Gamma_n(i,j) * norm(squeeze(En(i,j,:))) / norm(squeeze(CN(i,j,:)));
                end

                %% South
                if j > 1
                    aS(i,j) = -Gamma_s(i,j) * norm(squeeze(Es(i,j,:))) / norm(squeeze(CS(i,j,:)));
                end

                %% Central coefficient
                aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));

                %% Right-hand side
                bC(i,j) = S(i,j)*Vc(i,j) + bCorr(i,j);

            end
        end

    end
    function plotGeometry(X,Y,XC,YC,Se,Sw,Sn,Ss,Nx,Ny, Xe, Ye,Xn,Yn,CE)

        figure
        hold on
        axis equal
        box on

        %% ---- Plot grid lines ----
        [m,n] = size(X);

        for i = 1:m
            plot(X(i,:),Y(i,:),'k')
        end

        for j = 1:n
            plot(X(:,j),Y(:,j),'k')
        end

        %% ---- Plot nodes ----
        plot(X(:),Y(:),'r.','MarkerSize',10)

        %% ---- Plot centroids ----
        plot(XC(:),YC(:),'bo','MarkerSize',6,'LineWidth',1.5)

        %% ---- Plot surface vectors ----
        scale = 0.3;   % scale for visibility

        for i = 1:Nx
            for j = 1:Ny

                xc = XC(i,j);
                yc = YC(i,j);

                quiver(xc,yc,scale*Se(i,j,1),scale*Se(i,j,2),'r')
                quiver(xc,yc,scale*Sw(i,j,1),scale*Sw(i,j,2),'g')
                quiver(xc,yc,scale*Sn(i,j,1),scale*Sn(i,j,2),'b')
                quiver(xc,yc,scale*Ss(i,j,1),scale*Ss(i,j,2),'m')

            end
        end
        scale = 0.4;

        for i=1:Nx
            for j=1:Ny
                quiver(XC(i,j),YC(i,j),scale*CE(i,j,1),scale*CE(i,j,2),'k')
            end
        end
        plot(Xe(:),Ye(:),'gs','MarkerSize',6)
        plot(Xn(:),Yn(:),'ms','MarkerSize',6)
        title('Grid Geometry Check')
        xlabel('x')
        ylabel('y')

        legend('grid','nodes','centroids')

        hold off

    end
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
bc.east.type = 'robin';
bc.east.h = 15.0;
bc.east.Tinf = 300.0;

% choose for the south
bc.south.type = 'dirichlet';
bc.south.value = 320.0;

% choose for the north
bc.north.type = 'neumann';
bc.north.value = 0.0;

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

