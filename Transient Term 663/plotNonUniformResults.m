function plotNonUniformResults(XC, YC, T_be_save, T_cn_save, save_times, saveFolder)
%PLOTNONUNIFORMRESULTS Create clear comparison plots for BE and CN
%
% Inputs:
%   XC, YC      - cell centroid coordinates
%   T_be_save   - BE temperature field, size [Nx, Ny, Nt]
%   T_cn_save   - CN temperature field, size [Nx, Ny, Nt]
%   save_times  - times at which fields were saved, e.g. [10 30 60]
%   saveFolder  - folder where figures will be saved, e.g. 'figures'

    if nargin < 6 || isempty(saveFolder)
        saveFolder = 'figures';
    end

    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

    [Nx, Ny, Nt] = size(T_be_save);

    time_labels = arrayfun(@(t) sprintf('t = %ds', t), save_times, 'UniformOutput', false);

    % Centerlines
    mid_i = round(Nx/2);
    mid_j = round(Ny/2);

    y_line = YC(mid_i, :);   % vertical centerline coordinates
    x_line = XC(:, mid_j);   % horizontal centerline coordinates

    % Difference field
    diff_cn_be = T_cn_save - T_be_save;

    %% --------------------------------------------------------------------
    % 1) Centerline comparison plots
    % More telling than a single row because it shows both directions
    %% --------------------------------------------------------------------
    fig1 = figure('Name', 'NonUniform_Centerlines', 'Color', 'w');
    tl = tiledlayout(2, Nt, 'TileSpacing', 'compact', 'Padding', 'compact');

    for kk = 1:Nt
        % Vertical centerline
        nexttile
        plot(y_line, squeeze(T_be_save(mid_i,:,kk)), '-o', ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'BE');
        hold on
        plot(y_line, squeeze(T_cn_save(mid_i,:,kk)), '-s', ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'CN');
        grid on
        xlabel('y (m)')
        ylabel('T (K)')
        title(['Vertical centerline, ' time_labels{kk}])
        if kk == 1
            legend('Location', 'best')
        end

        % Horizontal centerline
        nexttile
        plot(x_line, squeeze(T_be_save(:,mid_j,kk)), '-o', ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'BE');
        hold on
        plot(x_line, squeeze(T_cn_save(:,mid_j,kk)), '-s', ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'CN');
        grid on
        xlabel('x (m)')
        ylabel('T (K)')
        title(['Horizontal centerline, ' time_labels{kk}])
    end

    title(tl, 'Temperature comparison along centerlines — non-uniform grid')
    exportgraphics(fig1, fullfile(saveFolder, 'nonuniform_centerline_comparison.png'), 'Resolution', 300);

    %% --------------------------------------------------------------------
    % 2) Contour plots with shared color limits
    % This makes BE and CN visually comparable
    %% --------------------------------------------------------------------
    for kk = 1:Nt
        fig2 = figure('Name', ['NonUniform_Contours_t' num2str(save_times(kk))], 'Color', 'w');
        tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        allT = [T_be_save(:,:,kk); T_cn_save(:,:,kk)];
        clim_range = [min(allT(:)), max(allT(:))];

        nexttile
        contourf(XC, YC, T_be_save(:,:,kk), 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(turbo)
        clim(clim_range)
        xlabel('x (m)')
        ylabel('y (m)')
        title(['BE, ' time_labels{kk}])

        nexttile
        contourf(XC, YC, T_cn_save(:,:,kk), 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(turbo)
        clim(clim_range)
        xlabel('x (m)')
        ylabel('y (m)')
        title(['CN, ' time_labels{kk}])

        title(tl, ['Non-uniform grid temperature contours — ' time_labels{kk}])

        exportgraphics(fig2, ...
            fullfile(saveFolder, sprintf('nonuniform_contours_t%d.png', save_times(kk))), ...
            'Resolution', 300);
    end

    %% --------------------------------------------------------------------
    % 3) CN - BE difference contours
    % Use symmetric color limits so positive/negative differences are fair
    %% --------------------------------------------------------------------
    for kk = 1:Nt
        fig3 = figure('Name', ['NonUniform_Difference_t' num2str(save_times(kk))], 'Color', 'w');

        current_diff = diff_cn_be(:,:,kk);
        diff_lim = max(abs(current_diff(:)));
        if diff_lim == 0
            diff_lim = 1;
        end

        contourf(XC, YC, current_diff, 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(parula)
        clim([-diff_lim, diff_lim])
        xlabel('x (m)')
        ylabel('y (m)')
        title(['CN - BE at ' time_labels{kk}])

        exportgraphics(fig3, ...
            fullfile(saveFolder, sprintf('nonuniform_diff_t%d.png', save_times(kk))), ...
            'Resolution', 300);
    end

    %% --------------------------------------------------------------------
    % 4) Difference along centerlines
    % This is more telling than just the contour because it shows where they differ
    %% --------------------------------------------------------------------
    fig4 = figure('Name', 'NonUniform_Centerline_Difference', 'Color', 'w');
    tl = tiledlayout(2, Nt, 'TileSpacing', 'compact', 'Padding', 'compact');

    for kk = 1:Nt
        nexttile
        plot(y_line, squeeze(diff_cn_be(mid_i,:,kk)), '-d', ...
            'LineWidth', 1.2, 'MarkerSize', 4);
        hold on
        yline(0, '--', 'LineWidth', 1.0);
        grid on
        xlabel('y (m)')
        ylabel('\DeltaT = T_{CN} - T_{BE} (K)')
        title(['Vertical centerline difference, ' time_labels{kk}])

        nexttile
        plot(x_line, squeeze(diff_cn_be(:,mid_j,kk)), '-d', ...
            'LineWidth', 1.2, 'MarkerSize', 4);
        hold on
        yline(0, '--', 'LineWidth', 1.0);
        grid on
        xlabel('x (m)')
        ylabel('\DeltaT = T_{CN} - T_{BE} (K)')
        title(['Horizontal centerline difference, ' time_labels{kk}])
    end

    title(tl, 'Centerline differences between CN and BE')
    exportgraphics(fig4, fullfile(saveFolder, 'nonuniform_centerline_difference.png'), 'Resolution', 300);

    %% --------------------------------------------------------------------
    % 5) Summary metrics vs time
    % Helpful if you want to discuss how far apart CN and BE are over time
    %% --------------------------------------------------------------------
    rms_diff = zeros(1, Nt);
    max_diff = zeros(1, Nt);
    mean_abs_diff = zeros(1, Nt);

    for kk = 1:Nt
        current_diff = diff_cn_be(:,:,kk);
        rms_diff(kk)      = sqrt(mean(current_diff.^2, 'all'));
        max_diff(kk)      = max(abs(current_diff), [], 'all');
        mean_abs_diff(kk) = mean(abs(current_diff), 'all');
    end

    fig5 = figure('Name', 'NonUniform_Difference_Summary', 'Color', 'w');
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    plot(save_times, rms_diff, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on
    xlabel('Time (s)')
    ylabel('RMS(CN - BE) (K)')
    title('RMS difference vs time')

    nexttile
    plot(save_times, max_diff, '-s', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on
    xlabel('Time (s)')
    ylabel('Max|CN - BE| (K)')
    title('Maximum absolute difference vs time')

    nexttile
    plot(save_times, mean_abs_diff, '-d', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on
    xlabel('Time (s)')
    ylabel('Mean|CN - BE| (K)')
    title('Mean absolute difference vs time')

    title(tl, 'Summary of differences between CN and BE')
    exportgraphics(fig5, fullfile(saveFolder, 'nonuniform_difference_summary.png'), 'Resolution', 300);

    %% --------------------------------------------------------------------
    % 6) Console summary
    %% --------------------------------------------------------------------
    fprintf('\n%-10s %-18s %-18s %-18s\n', ...
        'Time (s)', 'RMS(CN-BE)', 'Mean|CN-BE|', 'Max|CN-BE|');

    for kk = 1:Nt
        fprintf('%-10d %-18.6f %-18.6f %-18.6f\n', ...
            save_times(kk), rms_diff(kk), mean_abs_diff(kk), max_diff(kk));
    end
end