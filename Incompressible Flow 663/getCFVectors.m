%% the distance between the centers
function [CE,CW,CN,CS] = getCFVectors(XC,YC,Nx,Ny)

CE = zeros(Nx,Ny,2);
CW = zeros(Nx,Ny,2);
CN = zeros(Nx,Ny,2);
CS = zeros(Nx,Ny,2);

for i = 1:Nx
for j = 1:Ny

    %% ueast
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

    %% South
    if j > 1
        CS(i,j,1) = XC(i,j-1) -XC(i,j);
        CS(i,j,2) = YC(i,j-1) -YC(i,j);
    end

end
end

end
