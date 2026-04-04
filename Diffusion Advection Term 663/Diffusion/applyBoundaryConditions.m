function [aE,aW,aN,aS,aC,bC] = applyBoundaryConditions(aE,aW,aN,aS,aC,bC,Gamma,XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys, ...
    Se,Sw,Sn,Ss, Nx,Ny,bc) 
if nargin < 24
    bc.west.type = 'dirichlet';
    bc.west.value = 400.0;

    bc.east.type = 'robin';
    bc.east.h = 15.0;
    bc.east.Tinf = 300.0;

    bc.south.type = 'dirichlet';
    bc.south.value = 320.0;

    bc.north.type = 'neumann';
    bc.north.value = 0.0;
end

% West boundary.
for j = 1:Ny
    aW(1,j) = 0.0;
    area = norm(squeeze(Sw(1,j,:)));
    distance = hypot(Xw(1,j) - XC(1,j), Yw(1,j) - YC(1,j));

    switch lower(bc.west.type)
        case 'dirichlet'
            D = Gamma(1,j) * area / distance;
            aC(1,j) = aC(1,j) + D;
            bC(1,j) = bC(1,j) + D*bc.west.value;

        case 'neumann'
            bC(1,j) = bC(1,j) + Gamma(1,j)*area*bc.west.value;

        case 'robin'
            D = area / (distance/Gamma(1,j) + 1.0/bc.west.h);
            aC(1,j) = aC(1,j) + D;
            bC(1,j) = bC(1,j) + D*bc.west.Tinf;

        otherwise
            error('Please specify the BC.');
    end
end

% East boundary.
for j = 1:Ny
    aE(Nx,j) = 0.0;
    area = norm(squeeze(Se(Nx,j,:)));
    distance = hypot(Xe(Nx,j) - XC(Nx,j), Ye(Nx,j) - YC(Nx,j));

    switch lower(bc.east.type)
        case 'dirichlet'
            D = Gamma(Nx,j) * area / distance;
            aC(Nx,j) = aC(Nx,j) + D;
            bC(Nx,j) = bC(Nx,j) + D*bc.east.value;

        case 'neumann'
            bC(Nx,j) = bC(Nx,j) + Gamma(Nx,j)*area*bc.east.value;

        case 'robin'
            D = area / (distance/Gamma(Nx,j) + 1.0/bc.east.h);
            aC(Nx,j) = aC(Nx,j) + D;
            bC(Nx,j) = bC(Nx,j) + D*bc.east.Tinf;

        otherwise
            error('Please specify the BC.');
    end
end

% South boundary.
for i = 1:Nx
    aS(i,1) = 0.0;
    area = norm(squeeze(Ss(i,1,:)));
    distance = hypot(Xs(i,1) - XC(i,1), Ys(i,1) - YC(i,1));

    switch lower(bc.south.type)
        case 'dirichlet'
            D = Gamma(i,1) * area / distance;
            aC(i,1) = aC(i,1) + D;
            bC(i,1) = bC(i,1) + D*bc.south.value;

        case 'neumann'
            bC(i,1) = bC(i,1) + Gamma(i,1)*area*bc.south.value;

        case 'robin'
            D = area / (distance/Gamma(i,1) + 1.0/bc.south.h);
            aC(i,1) = aC(i,1) + D;
            bC(i,1) = bC(i,1) + D*bc.south.Tinf;

        otherwise
            error('Please specify the BC.');
    end
end

% North boundary.
for i = 1:Nx
    aN(i,Ny) = 0.0;

    area = norm(squeeze(Sn(i,Ny,:)));
    distance = hypot(Xn(i,Ny) - XC(i,Ny), Yn(i,Ny) - YC(i,Ny));

    switch lower(bc.north.type)
        case 'dirichlet'
            D = Gamma(i,Ny) * area / distance;
            aC(i,Ny) = aC(i,Ny) + D;
            bC(i,Ny) = bC(i,Ny) + D*bc.north.value;

        case 'neumann'
            bC(i,Ny) = bC(i,Ny) + Gamma(i,Ny)*area*bc.north.value;

        case 'robin'
            D = area / (distance/Gamma(i,Ny) + 1.0/bc.north.h);
            aC(i,Ny) = aC(i,Ny) + D;
            bC(i,Ny) = bC(i,Ny) + D*bc.north.Tinf;

        otherwise
            error('Please specify the BC.');
    end
end

end
