function [L2_error, Max_error, T_exact, errorField] = compareWithAnalytical(Tsol,X,Y,Nx,Ny,makePlots)

if nargin < 6
    makePlots = true;
end

xmin = min(X(:));
xmax = max(X(:));
Lx = xmax - xmin;

T_exact = zeros(Nx,Ny);
for j = 1:Ny
    for i = 1:Nx
        
        x_center = 0.25*( ...
            X(i,j) + X(i+1,j) + ...
            X(i,j+1) + X(i+1,j+1) );
        
        xhat = (x_center - xmin)/Lx;

        T_exact(i,j) = 400 - 100*xhat;

    end
end

%% error field
errorField = abs(Tsol - T_exact);

L2_error = sqrt(mean(errorField(:).^2));
Max_error = max(errorField(:));

fprintf('L2 Error = %e\n',L2_error)
fprintf('Max Error = %e\n',Max_error)

%% plots
if makePlots
    figure
    contourf(T_exact,20)
    colorbar
    title('Analytical Solution')

    figure
    contourf(Tsol,20)
    colorbar
    title('Numerical Solution')

    figure
    contourf(errorField,20)
    colorbar
    title('Error |T_{num}-T_{exact}|')
end

end
