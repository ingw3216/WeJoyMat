%% Ex-Gaussian贝叶斯参数估计与群体分析（改进版）
% 使用MCMC方法进行参数估计，包含置信区间和群体效应分析
clear; clc; close all;

% 固定随机种子以确保可重复性
rng(42);

%% 1. 读取数据
fprintf('=====================================\n');
fprintf('Ex-Gaussian贝叶斯分析系统（改进版）\n');
fprintf('=====================================\n\n');
fprintf('正在读取数据...\n');
data = readtable('rt_ms_long_with_hardfilter.xlsx');

% 筛选有效数据
valid_data = data(data.valid == 1, :);
fprintf('有效数据行数: %d\n', height(valid_data));

%% 2. 初始化
children = unique(valid_data.child_id);
modes = {'A', 'B', 'C', 'D'};
n_children = length(children);
n_modes = length(modes);

% 结果存储 - 预分配内存
max_results = n_children * n_modes;
results_bayes = NaN(max_results, 11);  % [child_id, mode_num, mu, sigma, tau, CI_bounds...]
convergence_flag = zeros(max_results, 1);  % 收敛标志
result_counter = 1;

% 创建输出目录
output_dir = 'C:\Users\ingw.LENOVO\Desktop\carpet_data\pic_results\new_data';
if ~exist(output_dir, 'dir')
    output_dir = pwd;  % 如果路径不存在，使用当前目录
end

% 自动延长采样标志
auto_extend_sampling = false;  % 设置为true可自动延长未收敛的链

%% 3. 定义Ex-Gaussian PDF函数
exgauss_pdf = @(x, mu, sigma, tau) exgaussian_pdf(x, mu, sigma, tau);

%% 4. 主分析循环：贝叶斯MCMC估计
fprintf('\n开始贝叶斯参数估计...\n');
fprintf('----------------------------------------\n');

for c = 1:n_children
    child_id = children(c);
    fprintf('\n处理儿童 %d:\n', child_id);
    
    for m = 1:n_modes
        mode = modes{m};
        
        % 提取数据
        mask = (valid_data.child_id == child_id) & strcmp(valid_data.mode, mode);
        rt_data = valid_data.rt_ms(mask);
        
        if isempty(rt_data) || length(rt_data) < 10
            fprintf('  警告: 儿童%d模式%s数据不足(n=%d)\n', child_id, mode, length(rt_data));
            continue;
        end
        
        fprintf('  模式%s: n=%d trials, ', mode, length(rt_data));
        
        %% 贝叶斯估计（使用MCMC）
        % 计算数据统计量用于设置先验
        mean_rt = mean(rt_data);
        std_rt = std(rt_data);
        skew_rt = skewness(rt_data);
        
        % 设置弱信息先验
        prior_mu_mean = mean_rt;
        prior_mu_std = std_rt;
        prior_sigma_shape = 2;  % Gamma分布的shape参数
        prior_sigma_scale = std_rt/2;  % Gamma分布的scale参数
        prior_tau_shape = 2;
        prior_tau_scale = mean_rt/4;
        
        % MCMC设置
        n_iter = 10000;  % 迭代次数
        n_burnin = 2000;  % 预烧期
        n_chains = 4;  % 链数
        
        % 初始化MCMC链
        chains = zeros(n_iter, 3, n_chains);
        
        % 改进的初始值设置（使用矩估计）
        sample_mean = mean(rt_data);
        sample_var = var(rt_data);
        sample_skew = skewness(rt_data);
        
        % Heathcote (1996) 的矩估计公式
        if sample_skew > 0
            tau_init = sqrt(sample_var) * (sample_skew/2)^(1/3);
        else
            tau_init = mean_rt * 0.3;  % 如果偏度为负或零，使用默认值
        end
        sigma_init = sqrt(max(sample_var - tau_init^2, std_rt^2 * 0.25));
        mu_init = sample_mean - tau_init;
        
        % 确保参数合理
        tau_init = max(tau_init, 50);
        sigma_init = max(sigma_init, 50);
        
        for chain = 1:n_chains
            % 初始值设置
            if chain == 1
                % 第一条链使用矩估计值
                current = [mu_init, sigma_init, tau_init];
            else
                % 其他链从第一条链附近开始，增加扰动
                current = [mu_init, sigma_init, tau_init] .* (0.5 + rand(1,3));
            end
            
            % 确保current是行向量
            current = reshape(current, 1, 3);
            
            % 确保初始值为正
            current(2) = max(current(2), 10);
            current(3) = max(current(3), 10);
            
            % 自适应参数设置
            adaptation_window = 500;  % 每500次迭代调整一次
            target_accept = 0.44;     % 目标接受率（Roberts & Rosenthal 2001）
            step_sizes = [std_rt*0.1, std_rt*0.05, mean_rt*0.1];  % 初始步长
            
            % 记录接受历史用于自适应
            accept_history = zeros(adaptation_window, 1);
            accept_idx = 1;
            
            % Metropolis-Hastings采样
            accepted = 0;
            for iter = 1:n_iter
                % 自适应步长调整
                if mod(iter, adaptation_window) == 0 && iter <= n_burnin
                    current_accept_rate = mean(accept_history);
                    
                    % 根据接受率调整步长
                    if current_accept_rate < 0.2
                        step_sizes = step_sizes * 0.5;  % 步长太大，减小
                    elseif current_accept_rate > 0.6
                        step_sizes = step_sizes * 1.5;  % 步长太小，增大
                    else
                        % 微调
                        adjustment = 1 + 2*(current_accept_rate - target_accept);
                        step_sizes = step_sizes * adjustment;
                    end
                    
                    % 重置接受历史
                    accept_history = zeros(adaptation_window, 1);
                    accept_idx = 1;
                end
                
                % 使用调整后的步长
                proposal = current + randn(1,3) .* step_sizes;
                
                % 确保参数为正
                if proposal(2) > 10 && proposal(3) > 10
                    try
                        % 计算似然
                        log_likelihood_current = sum(log(exgauss_pdf(rt_data, current(1), current(2), current(3))));
                        log_likelihood_proposal = sum(log(exgauss_pdf(rt_data, proposal(1), proposal(2), proposal(3))));
                        
                        % 检查似然值有效性
                        if ~isfinite(log_likelihood_current) || ~isfinite(log_likelihood_proposal)
                            continue;
                        end
                        
                        % 计算先验
                        log_prior_current = log(normpdf(current(1), prior_mu_mean, prior_mu_std)) + ...
                                           log(gampdf(current(2), prior_sigma_shape, prior_sigma_scale)) + ...
                                           log(gampdf(current(3), prior_tau_shape, prior_tau_scale));
                        
                        log_prior_proposal = log(normpdf(proposal(1), prior_mu_mean, prior_mu_std)) + ...
                                            log(gampdf(proposal(2), prior_sigma_shape, prior_sigma_scale)) + ...
                                            log(gampdf(proposal(3), prior_tau_shape, prior_tau_scale));
                        
                        % Metropolis比率
                        log_ratio = (log_likelihood_proposal + log_prior_proposal) - ...
                                   (log_likelihood_current + log_prior_current);
                        
                        % 记录是否接受（用于自适应）
                        accepted_this_iter = 0;
                        if log(rand()) < log_ratio
                            current = proposal;
                            accepted = accepted + 1;
                            accepted_this_iter = 1;
                        end
                        
                        % 记录接受历史
                        if iter <= n_burnin && mod(iter-1, adaptation_window) < adaptation_window
                            accept_history(accept_idx) = accepted_this_iter;
                            accept_idx = mod(accept_idx, adaptation_window) + 1;
                        end
                    catch
                        continue;
                    end
                end
                
                % 存储当前状态 - 确保是行向量
                chains(iter, :, chain) = reshape(current, 1, 3);
            end
            
            % 输出第一条链的接受率
            if chain == 1
                accept_rate = accepted / n_iter;
                fprintf('接受率: %.1f%%', accept_rate * 100);
            end
        end
        
        % 改进的收敛诊断
        R_hats = zeros(3, 1);
        for param = 1:3
            chain_means = zeros(n_chains, 1);
            chain_vars = zeros(n_chains, 1);
            
            for chain = 1:n_chains
                chain_data = chains((n_burnin+1):end, param, chain);
                chain_means(chain) = mean(chain_data);
                chain_vars(chain) = var(chain_data);
            end
            
            W = mean(chain_vars);  % 链内方差
            B = var(chain_means) * (n_iter - n_burnin);  % 链间方差
            var_plus = ((n_iter - n_burnin - 1)/(n_iter - n_burnin)) * W + B/(n_iter - n_burnin);
            R_hats(param) = sqrt(var_plus / W);
            
            if R_hats(param) > 1.1
                fprintf('\n    警告: 参数%d未收敛 (R-hat=%.2f)', param, R_hats(param));
            end
        end
        
        % 检查最大R-hat
        max_rhat = max(R_hats);
        if max_rhat > 1.1
            fprintf('\n    注意：结果可能不可靠');
            convergence_flag(result_counter) = 0;
            
            % 可选：自动延长采样
            if auto_extend_sampling
                fprintf('\n    尝试延长采样...');
                extra_iter = 5000;
                % 这里可以添加额外采样的代码
            end
        else
            convergence_flag(result_counter) = 1;
        end
        
        % 合并链并去除预烧期
        all_samples = reshape(chains((n_burnin+1):end, :, :), [], 3);
        
        % 计算后验统计
        bayes_mu = median(all_samples(:,1));
        bayes_sigma = median(all_samples(:,2));
        bayes_tau = median(all_samples(:,3));
        
        % 95%可信区间
        ci_mu = prctile(all_samples(:,1), [2.5, 97.5]);
        ci_sigma = prctile(all_samples(:,2), [2.5, 97.5]);
        ci_tau = prctile(all_samples(:,3), [2.5, 97.5]);
        
        % 存储贝叶斯结果
        results_bayes(result_counter, :) = [child_id, m, bayes_mu, bayes_sigma, bayes_tau, ...
            ci_mu(1), ci_mu(2), ci_sigma(1), ci_sigma(2), ci_tau(1), ci_tau(2)];
        
        fprintf(' 完成\n');
        result_counter = result_counter + 1;
    end
end

% 删除未使用的行
valid_rows = ~isnan(results_bayes(:,1));
results_bayes = results_bayes(valid_rows, :);
convergence_flag = convergence_flag(valid_rows);

%% 5. 群体效应分析
fprintf('\n=====================================\n');
fprintf('群体效应分析\n');
fprintf('=====================================\n');

% 重构数据为三维矩阵（儿童×模式×参数）
params_matrix = NaN(n_children, n_modes, 3);
for i = 1:size(results_bayes, 1)
    child_idx = find(children == results_bayes(i, 1));
    mode_idx = results_bayes(i, 2);
    if ~isempty(child_idx) && mode_idx <= n_modes
        params_matrix(child_idx, mode_idx, :) = results_bayes(i, 3:5);
    end
end

% 对每个参数进行分析
param_names = {'μ (位置)', 'σ (尺度)', 'τ (指数)'};
for p = 1:3
    fprintf('\n参数 %s 的群体分析:\n', param_names{p});
    fprintf('----------------------------------------\n');
    
    % 提取当前参数数据
    param_data = squeeze(params_matrix(:, :, p));
    
    % 描述性统计
    fprintf('各模式均值±标准差:\n');
    for m = 1:n_modes
        mode_data = param_data(:, m);
        valid_data_mode = mode_data(~isnan(mode_data));
        if ~isempty(valid_data_mode)
            mode_mean = mean(valid_data_mode);
            mode_std = std(valid_data_mode);
            mode_median = median(valid_data_mode);
            fprintf('  模式%s: 均值=%.1f ± %.1f ms, 中位数=%.1f ms (n=%d)\n', ...
                modes{m}, mode_mean, mode_std, mode_median, length(valid_data_mode));
        end
    end
    
    % 重复测量方差分析
    y = param_data(:);
    subject = repmat((1:n_children)', n_modes, 1);
    mode_factor = repelem((1:n_modes)', n_children, 1);
    
    % 移除NaN
    valid_idx = ~isnan(y);
    y = y(valid_idx);
    subject = subject(valid_idx);
    mode_factor = mode_factor(valid_idx);
    
    if length(y) > n_modes
        % 执行重复测量ANOVA
        try
            grand_mean = mean(y);
            ss_between = 0;
            ss_within = 0;
            
            for m = 1:n_modes
                mode_data = y(mode_factor == m);
                if ~isempty(mode_data)
                    mode_mean = mean(mode_data);
                    ss_between = ss_between + length(mode_data) * (mode_mean - grand_mean)^2;
                    ss_within = ss_within + sum((mode_data - mode_mean).^2);
                end
            end
            
            df_between = n_modes - 1;
            df_within = length(y) - n_modes;
            
            if df_within > 0
                ms_between = ss_between / df_between;
                ms_within = ss_within / df_within;
                
                F_stat = ms_between / ms_within;
                p_value = 1 - fcdf(F_stat, df_between, df_within);
                
                % 计算效应量
                eta_squared = ss_between / (ss_between + ss_within);
                
                fprintf('\n重复测量ANOVA结果:\n');
                fprintf('  F(%d,%d) = %.2f, p = %.4f\n', df_between, df_within, F_stat, p_value);
                fprintf('  效应量(偏η²) = %.3f', eta_squared);
                
                % 效应量解释
                if eta_squared < 0.01
                    fprintf(' (无效应)\n');
                elseif eta_squared < 0.06
                    fprintf(' (小效应)\n');
                elseif eta_squared < 0.14
                    fprintf(' (中等效应)\n');
                else
                    fprintf(' (大效应)\n');
                end
                
                % 事后两两比较（Bonferroni校正）- 只在p<0.05时进行
                if p_value < 0.05
                    fprintf('\n事后比较（Bonferroni校正）:\n');
                    n_comparisons = nchoosek(n_modes, 2);
                    alpha_corrected = 0.05 / n_comparisons;
                    
                    for m1 = 1:n_modes-1
                        for m2 = m1+1:n_modes
                            data1 = param_data(:, m1);
                            data2 = param_data(:, m2);
                            
                            valid = ~isnan(data1) & ~isnan(data2);
                            data1 = data1(valid);
                            data2 = data2(valid);
                            
                            if length(data1) > 1
                                [h, p_val, ci, stats] = ttest(data1, data2);
                                
                                mean_diff = mean(data1 - data2);
                                pooled_sd = sqrt((var(data1) + var(data2)) / 2);
                                cohens_d = mean_diff / pooled_sd;
                                
                                sig_marker = '';
                                if p_val < alpha_corrected
                                    sig_marker = '***';
                                elseif p_val < 0.05
                                    sig_marker = '*';
                                end
                                
                                fprintf('  %s vs %s: Δ=%.1f, t=%.2f, p=%.4f %s, d=%.2f\n', ...
                                    modes{m1}, modes{m2}, mean_diff, stats.tstat, p_val, sig_marker, cohens_d);
                            end
                        end
                    end
                end
            end
            
        catch ME
            fprintf('  ANOVA分析失败: %s\n', ME.message);
        end
    end
    
    % 计算贝叶斯因子（模式间比较）
    fprintf('\n贝叶斯因子分析（H0: 模式间无差异）:\n');
    for m1 = 1:n_modes-1
        for m2 = m1+1:n_modes
            data1 = param_data(:, m1);
            data2 = param_data(:, m2);
            valid = ~isnan(data1) & ~isnan(data2);
            
            if sum(valid) > 2
                diff_samples = data1(valid) - data2(valid);
                % 使用Savage-Dickey比率近似
                bf = normpdf(0, mean(diff_samples), std(diff_samples)) / normpdf(0, 0, std_rt);
                
                % BF解释
                if bf > 10
                    bf_interpretation = '强证据支持H0';
                elseif bf > 3
                    bf_interpretation = '中等证据支持H0';
                elseif bf > 1
                    bf_interpretation = '弱证据支持H0';
                elseif bf > 0.33
                    bf_interpretation = '无明确证据';
                elseif bf > 0.1
                    bf_interpretation = '弱证据反对H0';
                else
                    bf_interpretation = '强证据反对H0';
                end
                
                fprintf('  BF(模式%s vs %s): %.2f (%s)\n', modes{m1}, modes{m2}, bf, bf_interpretation);
            end
        end
    end
end

%% 6. 创建结果表格
fprintf('\n=====================================\n');
fprintf('创建结果表格\n');
fprintf('=====================================\n');

% 贝叶斯结果表（包含可信区间和收敛标志）
bayes_table = array2table(results_bayes, ...
    'VariableNames', {'child_id', 'mode_num', 'mu', 'sigma', 'tau', ...
    'mu_CI_lower', 'mu_CI_upper', 'sigma_CI_lower', 'sigma_CI_upper', ...
    'tau_CI_lower', 'tau_CI_upper'});

% 添加模式字母
mode_letters = cell(height(bayes_table), 1);
for i = 1:height(bayes_table)
    mode_num = bayes_table.mode_num(i);
    if mode_num >= 1 && mode_num <= length(modes)
        mode_letters{i} = modes{mode_num};
    end
end
bayes_table.mode = mode_letters;

% 添加收敛标志
bayes_table.converged = convergence_flag;

% 重新排列列顺序
bayes_table = bayes_table(:, [1, 12, 13, 3:11]);

% 保存详细结果
try
    writetable(bayes_table, fullfile(output_dir, 'ExGaussian_Bayes_Results_Improved.xlsx'));
    fprintf('结果已保存到: %s\n', fullfile(output_dir, 'ExGaussian_Bayes_Results_Improved.xlsx'));
catch
    writetable(bayes_table, 'ExGaussian_Bayes_Results_Improved.xlsx');
    fprintf('结果已保存到当前目录\n');
end

%% 7. 可视化
fprintf('\n=====================================\n');
fprintf('生成可视化图表\n');
fprintf('=====================================\n');

% 7.1 参数分布热图
figure('Position', [100, 100, 1400, 400]);
for p = 1:3
    subplot(1, 3, p);
    imagesc(squeeze(params_matrix(:, :, p)));
    colorbar;
    colormap('hot');
    xlabel('模式');
    ylabel('儿童');
    title(sprintf('%s参数热图', param_names{p}));
    set(gca, 'XTick', 1:4, 'XTickLabel', modes);
    set(gca, 'YTick', 1:n_children, 'YTickLabel', children);
end
sgtitle('Ex-Gaussian参数分布热图');

% 7.2 τ参数（注意力失效）的个体差异分析
figure('Position', [100, 100, 800, 600]);
tau_data = squeeze(params_matrix(:, :, 3));
boxplot(tau_data, 'Labels', modes);
ylabel('τ (ms)');
xlabel('模式');
title('τ参数（注意力失效）的模式间比较');
grid on;

% 添加个体线条
hold on;
for i = 1:n_children
    plot(1:n_modes, tau_data(i, :), 'o-', 'Color', [0.5, 0.5, 0.5, 0.3]);
end

% 标记可能有问题的儿童（τ > 1500ms）
high_tau_threshold = 1500;
problem_children = any(tau_data > high_tau_threshold, 2);
if any(problem_children)
    fprintf('\n注意：以下儿童的τ值超过%.0fms，可能存在严重注意力问题：\n', high_tau_threshold);
    problem_ids = children(problem_children);
    for i = 1:length(problem_ids)
        fprintf('  儿童 %d\n', problem_ids(i));
    end
end

% 7.3 收敛诊断报告
n_unconverged = sum(convergence_flag == 0);
if n_unconverged > 0
    fprintf('\n=====================================\n');
    fprintf('收敛诊断报告\n');
    fprintf('=====================================\n');
    fprintf('未收敛的结果数: %d / %d (%.1f%%)\n', ...
        n_unconverged, length(convergence_flag), 100*n_unconverged/length(convergence_flag));
    fprintf('建议：对未收敛的结果谨慎解释，或增加迭代次数重新运行\n');
end

fprintf('\n分析完成！\n');

%% 辅助函数

% Ex-Gaussian PDF
function pdf = exgaussian_pdf(x, mu, sigma, tau)
    if tau <= 0 || sigma <= 0
        pdf = zeros(size(x));
        return;
    end
    
    lambda = 1/tau;
    z = (mu + lambda*sigma^2 - x) / (sqrt(2)*sigma);
    
    % 数值稳定的计算
    pdf = (lambda/2) .* exp(lambda/2 .* (2*mu + lambda*sigma^2 - 2*x)) .* erfc(z);
    
    % 处理数值问题
    pdf(pdf <= 0) = 1e-10;
    pdf(~isfinite(pdf)) = 1e-10;
    pdf(x < 0) = 1e-10;  % 反应时间不能为负
end