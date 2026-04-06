function [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma, ...
    XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Se,Sw,Sn,Ss,Nx,Ny,bc,Fe,Fw,Fn,Fs)
%APPLYBOUNDARYCONDITIONS  Impose BCs on the assembled coefficient matrices.
%
%  Convection rule (upwind, outward-normal sign convention):
%    Inflow  face: F < 0  →  assembleConvectionUpwind put |F| into a_nb.
%                             applyBC zeros a_nb and moves |F|*T_bc into bC.
%                             aC is NOT touched (no double-counting).
%    Outflow face: F > 0  →  assembleConvectionUpwind already put F into aC.
%                             Nothing extra to do here for convection.

if nargin < 25
    Fe = []; Fw = []; Fn = []; Fs = [];
end

% ── WEST ─────────────────────────────────────────────────────────────────
for j = 1:Ny
    aW(1,j)  = 0;
    area     = norm(squeeze(Sw(1,j,:)));
    distance = hypot(Xw(1,j)-XC(1,j), Yw(1,j)-YC(1,j));
    Ff       = pickF(Fw, 1, j);

    switch lower(bc.west.type)
        case 'dirichlet'
            D        = Gamma(1,j)*area/distance;
            aC(1,j)  = aC(1,j) + D;
            bC(1,j)  = bC(1,j) + D*bc.west.value;
            if ~isempty(Ff) && Ff < 0          % inflow: move |Fw|*T_bc to RHS only
                bC(1,j) = bC(1,j) - Ff*bc.west.value;
            end

        case 'neumann'
            bC(1,j) = bC(1,j) + Gamma(1,j)*area*bc.west.value;

        case 'robin'
            D       = area / (distance/Gamma(1,j) + 1/bc.west.h);
            aC(1,j) = aC(1,j) + D;
            bC(1,j) = bC(1,j) + D*bc.west.Tinf;
            if ~isempty(Ff) && Ff < 0
                bC(1,j) = bC(1,j) - Ff*bc.west.Tinf;
            end

        case 'convective'
            if ~isempty(Ff) && Ff < 0
                bC(1,j) = bC(1,j) - Ff*bc.west.value;
            end

        otherwise
            error('Unknown BC type on west boundary: %s', bc.west.type);
    end
end

% ── EAST ─────────────────────────────────────────────────────────────────
for j = 1:Ny
    aE(Nx,j) = 0;
    area     = norm(squeeze(Se(Nx,j,:)));
    distance = hypot(Xe(Nx,j)-XC(Nx,j), Ye(Nx,j)-YC(Nx,j));
    Ff       = pickF(Fe, Nx, j);

    switch lower(bc.east.type)
        case 'dirichlet'
            D         = Gamma(Nx,j)*area/distance;
            aC(Nx,j)  = aC(Nx,j) + D;
            bC(Nx,j)  = bC(Nx,j) + D*bc.east.value;
            if ~isempty(Ff) && Ff < 0          % inflow from east (unusual but general)
                bC(Nx,j) = bC(Nx,j) - Ff*bc.east.value;
            end

        case 'neumann'
            bC(Nx,j) = bC(Nx,j) + Gamma(Nx,j)*area*bc.east.value;

        case 'robin'
            D         = area / (distance/Gamma(Nx,j) + 1/bc.east.h);
            aC(Nx,j)  = aC(Nx,j) + D;
            bC(Nx,j)  = bC(Nx,j) + D*bc.east.Tinf;
            if ~isempty(Ff) && Ff < 0
                bC(Nx,j) = bC(Nx,j) - Ff*bc.east.Tinf;
            end

        case 'convective'
            if ~isempty(Ff) && Ff < 0
                bC(Nx,j) = bC(Nx,j) - Ff*bc.east.value;
            end

        otherwise
            error('Unknown BC type on east boundary: %s', bc.east.type);
    end
end

% ── SOUTH ────────────────────────────────────────────────────────────────
for i = 1:Nx
    aS(i,1)  = 0;
    area     = norm(squeeze(Ss(i,1,:)));
    distance = hypot(Xs(i,1)-XC(i,1), Ys(i,1)-YC(i,1));
    Ff       = pickF(Fs, i, 1);

    switch lower(bc.south.type)
        case 'dirichlet'
            D       = Gamma(i,1)*area/distance;
            aC(i,1) = aC(i,1) + D;
            bC(i,1) = bC(i,1) + D*bc.south.value;
            if ~isempty(Ff) && Ff < 0
                bC(i,1) = bC(i,1) - Ff*bc.south.value;
            end

        case 'neumann'
            bC(i,1) = bC(i,1) + Gamma(i,1)*area*bc.south.value;

        case 'robin'
            D       = area / (distance/Gamma(i,1) + 1/bc.south.h);
            aC(i,1) = aC(i,1) + D;
            bC(i,1) = bC(i,1) + D*bc.south.Tinf;
            if ~isempty(Ff) && Ff < 0
                bC(i,1) = bC(i,1) - Ff*bc.south.Tinf;
            end

        case 'convective'
            if ~isempty(Ff) && Ff < 0
                bC(i,1) = bC(i,1) - Ff*bc.south.value;
            end

        otherwise
            error('Unknown BC type on south boundary: %s', bc.south.type);
    end
end

% ── NORTH ────────────────────────────────────────────────────────────────
for i = 1:Nx
    aN(i,Ny) = 0;
    area     = norm(squeeze(Sn(i,Ny,:)));
    distance = hypot(Xn(i,Ny)-XC(i,Ny), Yn(i,Ny)-YC(i,Ny));
    Ff       = pickF(Fn, i, Ny);

    switch lower(bc.north.type)
        case 'dirichlet'
            D        = Gamma(i,Ny)*area/distance;
            aC(i,Ny) = aC(i,Ny) + D;
            bC(i,Ny) = bC(i,Ny) + D*bc.north.value;
            if ~isempty(Ff) && Ff < 0
                bC(i,Ny) = bC(i,Ny) - Ff*bc.north.value;
            end

        case 'neumann'
            bC(i,Ny) = bC(i,Ny) + Gamma(i,Ny)*area*bc.north.value;

        case 'robin'
            D        = area / (distance/Gamma(i,Ny) + 1/bc.north.h);
            aC(i,Ny) = aC(i,Ny) + D;
            bC(i,Ny) = bC(i,Ny) + D*bc.north.Tinf;
            if ~isempty(Ff) && Ff < 0
                bC(i,Ny) = bC(i,Ny) - Ff*bc.north.Tinf;
            end

        case 'convective'
            if ~isempty(Ff) && Ff < 0
                bC(i,Ny) = bC(i,Ny) - Ff*bc.north.value;
            end

        otherwise
            error('Unknown BC type on north boundary: %s', bc.north.type);
    end
end

end

% ── helper ───────────────────────────────────────────────────────────────
function F = pickF(Farray, i, j)
    if isempty(Farray)
        F = [];
    else
        F = Farray(i,j);
    end
end