%% Find the properties at the centers, simple function
function [Gamma, S] = computeMaterialProperties(XC,YC,T,Nx,Ny)
Gamma = zeros(Nx,Ny);
S     = zeros(Nx,Ny);

%Loop over control volumes
for i = 1:Nx
for j = 1:Ny

    x = XC(i,j);
    y = YC(i,j);

    % Source term
    S(i,j) = T(i,j)*( 2*x - 0.2*y )/400;
    % Diffusion coefficients
    Gamma(i,j) = T(i,j)*( x^2 + exp(0.1*y) )/400;


end
end

end