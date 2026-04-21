function plotTransientResults(XC, YC, T_be_save, T_cn_save, T_anal_save, save_times)
%PLOTTRANSIENTRESULTS Create comparison, contour, and error plots
%
% Inputs:
%   XC, YC        - cell centroid coordinates
%   T_be_save     - backward Euler temperatures, size [Nx Ny Nt]
%   T_cn_save     - Crank-Nicolson temperatures, size [Nx Ny Nt]
%   T_anal_save   - analytical temperatures, size [Nx Ny Nt]
%   save_times    - vector of saved times, e.g. [10 30 60]

    time_labels = arrayfun(@(t) sprintf('t = %ds', t), save_times, 'UniformOutput', false);

    [Nx, Ny, Nt] = size(T_be_save);
    mid_i = round(Nx/2);
    mid_j = round(Ny/2);

    x_line = XC(:, mid_j);
    y_line = YC(mid_i, :);

    % Errors
    err_be = T_be_save - T_anal_save;
    err_cn = T_cn_save - T_anal_save;
    abs_err_be = abs(err_be);
    abs_err_cn = abs(err_cn);

    %% 1) Centerline comparison plots
    fig1 = figure('Name', 'Centerline Comparisons', 'Color', 'w');
    tl = tiledlayout(2, Nt, 'TileSpacing', 'compact', 'Padding', 'compact');

    for kk = 1:Nt
        % Vertical centerline
        nexttile
        plot(y_line, squeeze(T_be_save(mid_i,:,kk)),  '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'BE');
        hold on
        plot(y_line, squeeze(T_cn_save(mid_i,:,kk)),  '-s', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'CN');
        plot(y_line, squeeze(T_anal_save(mid_i,:,kk)), '--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
        grid on
        xlabel('y (m)')
        ylabel('T (K)')
        title(['Vertical centerline, ' time_labels{kk}])

        if kk == 1
            legend('Location', 'best')
        end

        % Horizontal centerline
        nexttile
        plot(x_line, squeeze(T_be_save(:,mid_j,kk)),  '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'BE');
        hold on
        plot(x_line, squeeze(T_cn_save(:,mid_j,kk)),  '-s', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'CN');
        plot(x_line, squeeze(T_anal_save(:,mid_j,kk)), '--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
        grid on
        xlabel('x (m)')
        ylabel('T (K)')
        title(['Horizontal centerline, ' time_labels{kk}])
    end

    title(tl, 'Temperature comparison along centerlines')
    exportgraphics(fig1, 'centerline_comparison.png', 'Resolution', 300);

    %% 2) Temperature contour comparisons
    for kk = 1:Nt
        fig2 = figure('Name', ['Contours_' num2str(save_times(kk))], 'Color', 'w');
        tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        allT = cat(1, T_be_save(:,:,kk), T_cn_save(:,:,kk), T_anal_save(:,:,kk));
        clim_range = [min(allT(:)), max(allT(:))];

        nexttile
        contourf(XC, YC, T_be_save(:,:,kk), 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(turbo)
        clim(clim_range)
        xlabel('x (m)'); ylabel('y (m)')
        title(['BE, ' time_labels{kk}])

        nexttile
        contourf(XC, YC, T_cn_save(:,:,kk), 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(turbo)
        clim(clim_range)
        xlabel('x (m)'); ylabel('y (m)')
        title(['CN, ' time_labels{kk}])

        nexttile
        contourf(XC, YC, T_anal_save(:,:,kk), 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(turbo)
        clim(clim_range)
        xlabel('x (m)'); ylabel('y (m)')
        title(['Analytical, ' time_labels{kk}])

        title(tl, ['Temperature contours at ' time_labels{kk}])
        exportgraphics(fig2, sprintf('temperature_contours_t%d.png', save_times(kk)), 'Resolution', 300);
    end

    %% 3) Absolute error contour comparisons
    for kk = 1:Nt
        fig3 = figure('Name', ['AbsError_' num2str(save_times(kk))], 'Color', 'w');
        tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        errmax = max([abs_err_be(:,:,kk), abs_err_cn(:,:,kk)], [], 'all');
        clim_range = [0, errmax];

        nexttile
        contourf(XC, YC, abs_err_be(:,:,kk), 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(parula)
        clim(clim_range)
        xlabel('x (m)'); ylabel('y (m)')
        title(['|BE - Analytical|, ' time_labels{kk}])

        nexttile
        contourf(XC, YC, abs_err_cn(:,:,kk), 24, 'LineColor', 'none');
        axis equal tight
        colorbar
        colormap(parula)
        clim(clim_range)
        xlabel('x (m)'); ylabel('y (m)')
        title(['|CN - Analytical|, ' time_labels{kk}])

        title(tl, ['Absolute error contours at ' time_labels{kk}])
        exportgraphics(fig3, sprintf('absolute_error_contours_t%d.png', save_times(kk)), 'Resolution', 300);
    end

    %% 4) Error along centerlines
    fig4 = figure('Name', 'Centerline Errors', 'Color', 'w');
    tl = tiledlayout(2, Nt, 'TileSpacing', 'compact', 'Padding', 'compact');

    for kk = 1:Nt
        nexttile
        plot(y_line, squeeze(abs_err_be(mid_i,:,kk)), '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', '|BE - Analytical|');
        hold on
        plot(y_line, squeeze(abs_err_cn(mid_i,:,kk)), '-s', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', '|CN - Analytical|');
        grid on
        xlabel('y (m)')
        ylabel('|Error| (K)')
        title(['Vertical centerline error, ' time_labels{kk}])
        if kk == 1
            legend('Location', 'best')
        end

        nexttile
        plot(x_line, squeeze(abs_err_be(:,mid_j,kk)), '-o', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', '|BE - Analytical|');
        hold on
        plot(x_line, squeeze(abs_err_cn(:,mid_j,kk)), '-s', 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', '|CN - Analytical|');
        grid on
        xlabel('x (m)')
        ylabel('|Error| (K)')
        title(['Horizontal centerline error, ' time_labels{kk}])
    end

    title(tl, 'Absolute error along centerlines')
    exportgraphics(fig4, 'centerline_absolute_error.png', 'Resolution', 300);

    %% 5) Summary metrics vs time
    rms_be = zeros(1, Nt);
    rms_cn = zeros(1, Nt);
    max_be = zeros(1, Nt);
    max_cn = zeros(1, Nt);

    for kk = 1:Nt
        rms_be(kk) = sqrt(mean(err_be(:,:,kk).^2, 'all'));
        rms_cn(kk) = sqrt(mean(err_cn(:,:,kk).^2, 'all'));
        max_be(kk) = max(abs(err_be(:,:,kk)), [], 'all');
        max_cn(kk) = max(abs(err_cn(:,:,kk)), [], 'all');
    end

    fig5 = figure('Name', 'Error Summary', 'Color', 'w');
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    plot(save_times, rms_be, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'BE');
    hold on
    plot(save_times, rms_cn, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'CN');
    grid on
    xlabel('Time (s)')
    ylabel('RMS error (K)')
    title('RMS error vs time')
    legend('Location', 'best')

    nexttile
    plot(save_times, max_be, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'BE');
    hold on
    plot(save_times, max_cn, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'CN');
    grid on
    xlabel('Time (s)')
    ylabel('Max absolute error (K)')
    title('Max absolute error vs time')
    legend('Location', 'best')

    title(tl, 'Error summary metrics')
    exportgraphics(fig5, 'error_summary.png', 'Resolution', 300);

    %% 6) Console summary table
    fprintf('\n%-10s %-15s %-15s %-15s %-15s\n', ...
        'Time (s)', 'RMS BE', 'RMS CN', 'Max |BE|', 'Max |CN|');

    for kk = 1:Nt
        fprintf('%-10d %-15.6f %-15.6f %-15.6f %-15.6f\n', ...
            save_times(kk), rms_be(kk), rms_cn(kk), max_be(kk), max_cn(kk));
    end
end