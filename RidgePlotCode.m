%% Multisensory Interactive Carpet System for Children with Autism - Reaction Time Ridge Plot
% Function: Generate ridge plots for all 11 children's reaction time distributions across four sensory modes

clear; clc; close all;

%% 1. Data Reading and Preprocessing
% Read reaction time data
rt_filename = 'rt_ms_long_with_hardfilter.xlsx';
rt_data = readtable(rt_filename);

% Read Ex-Gaussian parameters
exgauss_filename = 'ExGaussian_Bayes_Results_Improved.xlsx';
exgauss_data = readtable(exgauss_filename);

% Filter valid data (valid=1)
valid_data = rt_data(rt_data.valid == 1, :);

% Define sensory modes and colors
modes = {'A', 'B', 'C', 'D'};
mode_names = {'Visual Only', 'Visual+Auditory', 'Visual+Tactile', 'Visual+Auditory+Tactile'};
mode_colors = [0.2 0.4 0.8;   % A: Blue
               0.8 0.2 0.2;   % B: Red
               0.2 0.7 0.3;   % C: Green
               0.9 0.6 0.1];  % D: Orange

% Create save directory
save_dir = 'C:\Users\ingw.LENOVO\Desktop\carpet_data\pic_results\new_data';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% 2. Loop over all 11 children
children = 1:11;  % 如果想根据数据自动识别，可用：children = unique(valid_data.child_id)';

% Common params
paramLabelSize    = 11;   % μ / μ+τ 字号
paramLineWidth    = 3;    % 竖线线宽
paramLabelYOffset = 0.1;  % μ/μ+τ 文字上移

x_min = 0;
x_max = 9000;
x_range = linspace(x_min, x_max, 200);

vertical_spacing = 1.0;
bottom_offset = -0.5;  % 抵消系统自动的空间

for child = children
    % Extract data for current child
    current_data = valid_data(valid_data.child_id == child, :);
    child_label = sprintf('C%d', child);

    % Prepare data containers
    ridge_data = cell(4, 1);
    ridge_labels = cell(4, 1);
    exgauss_params = struct('mu', zeros(4,1), 'sigma', zeros(4,1), 'tau', zeros(4,1));

    for mode_idx = 1:4
        mode = modes{mode_idx};
        mode_data = current_data(strcmp(current_data.mode, mode), :);

        % Store reaction time data
        ridge_data{mode_idx} = mode_data.rt_ms;
        ridge_labels{mode_idx} = sprintf('%s-Mode %s', child_label, mode);

        % Get Ex-Gaussian parameters
        exgauss_row = exgauss_data(exgauss_data.child_id == child & strcmp(exgauss_data.mode, mode), :);
        if ~isempty(exgauss_row)
            exgauss_params.mu(mode_idx)    = exgauss_row.mu;
            exgauss_params.sigma(mode_idx) = exgauss_row.sigma;
            exgauss_params.tau(mode_idx)   = exgauss_row.tau;
        end
    end

    %% 3. Create Ridge Plot (new figure for each child)
    fig = figure('Units','pixels','Position',[100, 100, 1200, 750]);
    hold on;

    for i = 1:4
        if ~isempty(ridge_data{i})
            % Kernel density estimation
            [density, x_density] = ksdensity(ridge_data{i}, x_range);

            % Normalize density values
            density = density / max(density) * 0.8;

            % Vertical position
            y_offset = (4 - i + 1) * vertical_spacing + bottom_offset;

            % Only draw from 200ms onwards
            valid_idx = x_density >= 200;
            x_plot = x_density(valid_idx);
            density_plot = density(valid_idx);

            % Filled area
            fill([x_plot, fliplr(x_plot)], ...
                 [y_offset + density_plot, y_offset * ones(1, length(density_plot))], ...
                 mode_colors(i, :), 'FaceAlpha', 0.6, 'EdgeColor', mode_colors(i, :), ...
                 'LineWidth', 1.5);

            % Baseline
            plot(x_range, y_offset * ones(size(x_range)), 'k-', 'LineWidth', 0.5);

            % Left label
            text(-200, y_offset + 0.2, ridge_labels{i}, ...
                 'FontSize', 11, 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'middle', 'FontWeight', 'bold');

            % Stats on the right
            n_points = length(ridge_data{i});
            if exgauss_params.mu(i) > 0
                mu_plus_tau = exgauss_params.mu(i) + exgauss_params.tau(i);
                stats_text = sprintf('μ = %.0f ms\nμ + τ = %.0f ms\nn = %d', ...
                    exgauss_params.mu(i), mu_plus_tau, n_points);
            else
                stats_text = sprintf('n = %d', n_points);
            end

            text(x_max + 200, y_offset + 0.2, stats_text, ...
                 'FontSize', 8, 'HorizontalAlignment', 'left', ...
                 'VerticalAlignment', 'middle', 'Color', [0.4 0.4 0.4]);

            % Ex-Gaussian parameter lines
            if exgauss_params.mu(i) > 0
                mu_density = interp1(x_density, density, exgauss_params.mu(i), 'linear', 0);
                exgauss_mean = exgauss_params.mu(i) + exgauss_params.tau(i);
                mean_density = interp1(x_density, density, exgauss_mean, 'linear', 0);

                % mu line
                plot([exgauss_params.mu(i), exgauss_params.mu(i)], ...
                     [y_offset, y_offset + mu_density], '--', 'Color', 'k', ...
                     'LineWidth', paramLineWidth);

                % mu + tau line
                plot([exgauss_mean, exgauss_mean], ...
                     [y_offset, y_offset + mean_density], '-', 'Color', 'k', ...
                     'LineWidth', paramLineWidth);

                % labels
                text(exgauss_params.mu(i),  y_offset + mu_density  + paramLabelYOffset, 'μ', ...
                    'FontSize', paramLabelSize, 'HorizontalAlignment','center', ...
                    'Color','k', 'FontWeight','bold');

                text(exgauss_mean,          y_offset + mean_density + paramLabelYOffset, 'μ+τ', ...
                    'FontSize', paramLabelSize, 'HorizontalAlignment','center', ...
                    'Color','k', 'FontWeight','bold');
            end
        end
    end

    %% 4. Chart Formatting - 混合方案：保留刻度功能，叠加箭头
    xlim([-300, x_max + 500]);
    ylim([0, 5]);

    ax = gca;
    ax.Box = 'off';
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.LineWidth = 1.5;
    ax.XColor = [0.5 0.5 0.5];
    ax.YColor = [0.5 0.5 0.5];
    ax.XTick = 0:1000:9000;
    ax.YTick = [];

    % annotation 需要使用归一化坐标
    ax_pos = ax.Position;
    x_lims = xlim; y_lims = ylim;

    x0_norm = ax_pos(1) + ax_pos(3) * (0 - x_lims(1)) / (x_lims(2) - x_lims(1));
    y0_norm = ax_pos(2) + ax_pos(4) * (0 - y_lims(1)) / (y_lims(2) - y_lims(1));
    x_end_norm = ax_pos(1) + ax_pos(3) * 0.97;
    y_end_norm = ax_pos(2) + ax_pos(4) * 0.95;

    annotation('arrow', [x0_norm, x_end_norm], [y0_norm, y0_norm], ...
               'HeadStyle', 'vback2', 'HeadLength', 8, 'HeadWidth', 8, ...
               'LineWidth', 1.5, 'Color', 'k');

    annotation('arrow', [x0_norm, x0_norm], [y0_norm, y_end_norm], ...
               'HeadStyle', 'vback2', 'HeadLength', 8, 'HeadWidth', 8, ...
               'LineWidth', 1.5, 'Color', 'k');

    xlabel('Reaction Time (ms)', 'FontSize', 15, 'FontWeight', 'bold');

    % 标题随儿童编号自动更新
    title(sprintf('Reaction Time Ridge Plot for %s', child_label), 'FontSize', 18, 'FontWeight', 'bold');

    set(gca, 'Color', [0.98 0.98 0.98]);

    % 图例
    legend_handles = gobjects(1,4);
    for i = 1:4
        legend_handles(i) = patch(NaN, NaN, mode_colors(i, :), 'FaceAlpha', 0.6);
    end
    legend(legend_handles, {sprintf('Mode A: %s', mode_names{1}), ...
                           sprintf('Mode B: %s', mode_names{2}), ...
                           sprintf('Mode C: %s', mode_names{3}), ...
                           sprintf('Mode D: %s', mode_names{4})}, ...
           'Location', 'northeast', 'FontSize', 12);

    %% 5. Save Figure (name per child)
    filename = sprintf('Ridge_Plot_%s.png', child_label);
    filepath = fullfile(save_dir, filename);

    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, filepath, '-dpng', '-r300');

    fprintf('Ridge plot for %s saved: %s\n', child_label, filepath);
    hold off; close(fig);
end
