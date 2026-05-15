%% same as notes
function [dTdx,dTdy] = computeCellGradient(Te,Tw,Tn,Ts,Se,Sw,Sn,Ss,Vc,Nx,Ny)
dTdx = zeros(Nx,Ny);
dTdy = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny

    dTdx(i,j) = (Te(i,j)*Se(i,j,1)+Tw(i,j)*Sw(i,j,1)+Tn(i,j)*Sn(i,j,1)+Ts(i,j)*Ss(i,j,1))/Vc(i,j);
    dTdy(i,j) = (Te(i,j)*Se(i,j,2)+Tw(i,j)*Sw(i,j,2)+Tn(i,j)*Sn(i,j,2)+Ts(i,j)*Ss(i,j,2))/Vc(i,j);

end
end

end