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