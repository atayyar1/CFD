%% Face Centers, taken to have 1 center
function [Xe,Ye,Xw,Yw,Xn,Yn,Xs,Ys] = getFaceCenters(X,Y,Nx,Ny)

Xe = zeros(Nx,Ny); Ye = zeros(Nx,Ny);

Xw = zeros(Nx,Ny); Yw = zeros(Nx,Ny);

Xn = zeros(Nx,Ny); Yn = zeros(Nx,Ny);

Xs = zeros(Nx,Ny); Ys = zeros(Nx,Ny);

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