%% please note that this was %fully% generated with AI tool
function plotFlowResults(XC, YC, u, v, p, T, Nx, Ny)
% PLOTFLOWRESULTS  Produce and save four plots for one grid resolution.
%
% Saves:
%   figures/figures_NxXNy.pdf        -> all plots in one PDF file
%   figures/velocity_NxXNy.fig       -> velocity figure
%   figures/streamlines_NxXNy.fig    -> streamlines figure
%   figures/temperature_NxXNy.fig    -> temperature figure
%   figures/pressure_NxXNy.fig       -> pressure figure
%
% XC, YC, u, v, p, T are all [Nx, Ny].
% MATLAB contourf / streamslice / streamline expect [Ny, Nx],
% so fields are transposed before plotting.

label = sprintf('%d\\times%d', Nx, Ny);
label_file = sprintf('%dx%d', Nx, Ny);

% =========================================================================
% Create output folder safely
% =========================================================================
outDir = fullfile(pwd, 'figures');

if ~isfolder(outDir)
    mkdir(outDir);
end

pdfFile = fullfile(outDir, ['figures_' label_file '.pdf']);

% Delete old PDF for this grid resolution if it already exists
if exist(pdfFile, 'file')
    delete(pdfFile);
end

% =========================================================================
% Transpose data to [Ny, Nx] for MATLAB plotting
% =========================================================================
x  = XC(:,1)';        % [1, Nx] x-coordinates of cell centres
y  = YC(1,:);         % [1, Ny] y-coordinates of cell centres

U  = u';              % [Ny, Nx]
V  = v';
P  = p';
TT = T';

speed = sqrt(U.^2 + V.^2);

% =========================================================================
% Figure 1 — Velocity magnitude + quiver
% =========================================================================
fig1 = figure('Name', ['Velocity  ' label], 'NumberTitle', 'off');

contourf(x, y, speed, 24, 'LineColor', 'none');
colorbar;
colormap(gca, 'turbo');

title(['Velocity magnitude [m/s]  —  ' label]);
xlabel('x [m]');
ylabel('y [m]');
axis equal tight;

% Skip arrows on fine grids by subsampling
step = max(1, round(Nx/20));

hold on
quiver(x(1:step:end), y(1:step:end), ...
       U(1:step:end, 1:step:end), ...
       V(1:step:end, 1:step:end), ...
       0.8, 'w', 'LineWidth', 0.8);
hold off

savefig(fig1, fullfile(outDir, ['velocity_' label_file '.fig']));
exportgraphics(fig1, pdfFile, 'Append', true);
% =========================================================================
% Figure 2 — Velocity vectors / quiver
% =========================================================================
fig2 = figure('Name', ['Velocity vectors  ' label], 'NumberTitle', 'off');

contourf(x, y, speed, 40, 'LineColor', 'none');
colormap(gca, 'hot');
colorbar;
hold on;

% Subsample vectors so the plot is readable
stepX = max(1, round(Nx/25));
stepY = max(1, round(Ny/12));

quiver(x(1:stepX:end), y(1:stepY:end), ...
       U(1:stepY:end, 1:stepX:end), ...
       V(1:stepY:end, 1:stepX:end), ...
       1.2, 'w', 'LineWidth', 0.8);

title(['Velocity vectors  —  ' label]);
xlabel('x [m]');
ylabel('y [m]');
axis equal tight;
hold off

savefig(fig2, fullfile(outDir, ['quiver_' label_file '.fig']));
exportgraphics(fig2, pdfFile, 'Append', true);

% =========================================================================
% Figure 3 — Temperature
% =========================================================================
fig3 = figure('Name', ['Temperature  ' label], 'NumberTitle', 'off');

contourf(x, y, TT, 24, 'LineColor', 'none');
colorbar;
colormap(gca, 'hot');

title(['Temperature [K]  —  ' label]);
xlabel('x [m]');
ylabel('y [m]');
axis equal tight;

% Isotherms overlay
hold on
contour(x, y, TT, 12, 'k', 'LineWidth', 0.5);
hold off

savefig(fig3, fullfile(outDir, ['temperature_' label_file '.fig']));
exportgraphics(fig3, pdfFile, 'Append', true);

% =========================================================================
% Figure 4 — Pressure gauge, relative to min
% =========================================================================
fig4 = figure('Name', ['Pressure  ' label], 'NumberTitle', 'off');

P_gauge = P - min(P(:));

contourf(x, y, P_gauge, 24, 'LineColor', 'none');
colorbar;
colormap(gca, 'parula');

title(['Pressure (gauge) [Pa]  —  ' label]);
xlabel('x [m]');
ylabel('y [m]');
axis equal tight;

savefig(fig4, fullfile(outDir, ['pressure_' label_file '.fig']));
exportgraphics(fig4, pdfFile, 'Append', true);

% =========================================================================
% Done
% =========================================================================
fprintf('\nFigures saved successfully.\n');
fprintf('Folder: %s\n', outDir);
fprintf('PDF:    %s\n\n', pdfFile);

end