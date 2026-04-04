function [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T,Gamma,bc,XC,YC, ...
    Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny)
%APPLY BOUNDARY FACE TEMPERATURES
%For gradient evaluation.

for j = 1:Ny
    distance = hypot(Xw(1,j) - XC(1,j), Yw(1,j) - YC(1,j));
    Tw(1,j) = boundaryFaceValue(T(1,j),Gamma(1,j),distance,bc.west);
end

for j = 1:Ny
    distance = hypot(Xe(Nx,j) - XC(Nx,j), Ye(Nx,j) - YC(Nx,j));
    Te(Nx,j) = boundaryFaceValue( ...
        T(Nx,j),Gamma(Nx,j),distance,bc.east);
end

for i = 1:Nx
    distance = hypot(Xs(i,1) - XC(i,1), Ys(i,1) - YC(i,1));
    Ts(i,1) = boundaryFaceValue( ...
        T(i,1),Gamma(i,1),distance,bc.south);
end

for i = 1:Nx
    distance = hypot(Xn(i,Ny) - XC(i,Ny), Yn(i,Ny) - YC(i,Ny));
    Tn(i,Ny) = boundaryFaceValue( ...
        T(i,Ny),Gamma(i,Ny),distance,bc.north);
end

end

function Tf = boundaryFaceValue(Tp,GammaP,distance,bcSide)

switch lower(bcSide.type)
    case 'dirichlet'
        Tf = bcSide.value;

    case 'neumann'
        Tf = Tp + distance*bcSide.value;

    case 'robin'
        Tf = (GammaP*Tp/distance + bcSide.h*bcSide.Tinf) / ...
             (GammaP/distance + bcSide.h);
    case 'convective'
        % inflow: use BC value, outflow: use cell value (zero gradient)
        if isfield(bcSide, 'value')
            Tf = bcSide.value;
        else
            Tf = Tp;  % outflow boundary
        end
    otherwise
        error('Unknown BC type. Use dirichlet, neumann, or robin.');
end

end
