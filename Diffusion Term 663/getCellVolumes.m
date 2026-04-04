%% Second function, Cell volumes, Simple as well
function elemVolume = getCellVolumes(X,Y)

[m,n] = size(X); Nx = m-1; Ny = n-1; elemVolume = zeros(Nx,Ny);
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

        % Cell area, the trapizoid trick
        elemVolume(i,j) = A1 + A2;

    end
end

end