
function [Se,Sw,Sn,Ss] = getSurfaceVectors(X,Y,Nx,Ny)

Se = zeros(Nx,Ny,2);
Sw = zeros(Nx,Ny,2);
Sn = zeros(Nx,Ny,2);
Ss = zeros(Nx,Ny,2);

for i = 1:Nx
    for j = 1:Ny

        %% East Face
        dx = X(i+1,j+1) - X(i+1,j); dy = Y(i+1,j+1) - Y(i+1,j);

        Se(i,j,1) = dy;
        Se(i,j,2) = -dx;

        %% west face
        dx = X(i,j+1)- X(i,j);
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
%% First I had them all positive, then changed it to face the other direction