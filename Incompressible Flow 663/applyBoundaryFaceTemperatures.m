function [Te,Tw,Tn,Ts] = applyBoundaryFaceTemperatures(Te,Tw,Tn,Ts,T,Gamma,bc, ...
    XC,YC,Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys,Nx,Ny,Fe,Fw,Fn,Fs)
%APPLYBOUNDARYFACETEMPERATURES  I Set face temperatures on boundary faces.
%  Used for gradient reconstruction, non-orthogonal correction, and I changed a bit in it for the SMART.
%
%  Fe,Fw,Fn,Fs are optional.  They are only needed when any BC type is
%  'convective', but I still added it
if nargin < 23 %I added this just in case
    Fe = []; Fw = []; Fn = []; Fs = [];
end
 for j = 1: Ny
   distance = hypot(Xw(1,j)-XC(1,j), Yw(1,j)-YC(1,j));
    Tw(1,j) = bfaceValue(T(1,j), Gamma(1,j), distance, bc.west, pickF(Fw,1,j));
end
 for j = 1: Ny
    distance = hypot(Xe(Nx,j)-XC(Nx,j), Ye(Nx,j)-YC(Nx,j));
    Te(Nx,j) = bfaceValue(T(Nx,j), Gamma(Nx,j), distance, bc.east, pickF(Fe,Nx,j));
end
 for i = 1: Nx 
    distance = hypot(Xs(i,1)-XC(i,1), Ys(i,1)-YC(i,1));   Ts(i,1) = bfaceValue(T(i,1), Gamma(i,1), distance, bc.south, pickF(Fs,i,1));
end
 for i = 1: Nx
   distance = hypot(Xn(i,Ny)-XC(i,Ny), Yn(i,Ny)-YC(i,Ny));
    Tn(i,Ny) = bfaceValue(T(i,Ny), Gamma(i,Ny), distance, bc.north, pickF(Fn,i,Ny));
end
end
 function Tf = bfaceValue(Tp, GammaP, distance, bcSide, Ff)
 switch lower(bcSide.type)
    case 'dirichlet'
        Tf = bcSide.value;

    case 'neumann'
        % bcSide.value is the outward-normal gradient dT/dn
        Tf = Tp + distance * bcSide.value;
    case 'robin'
        % harmonic blend of cell and ambient
        Tf = (GammaP*Tp/distance + bcSide.h*bcSide.Tinf) / ...
            (GammaP/distance + bcSide.h);

    case 'convective' 
        if ~isempty(Ff) && Ff < 0   
              Tf = bcSide.value;
       else                         
       Tf = Tp;
        end

    otherwise
        error('Unknown BC type: %s', bcSide.type);
end
end


% helper
function F = pickF(Farray, i, j)
if isempty(Farray)
    F = [];
else
    F = Farray(i,j);
end
end