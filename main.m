%
% ------------------------------------------------------------------------%
% Refer to:
% Yang C, Yang W, Qiu X, et al. Cognitive Radar Waveform Design Method
% under the Joint Constraints of Transmit Energy and Spectrum Bandwidth
% [J]. Remote Sensing, 2023, 15(21): 5187.
% ------------------------------------------------------------------------%
clc;clear;
close all;
rng default

bandwidth = 10e3; % 10k
pulse_width = 1e-3; % 1 ms
Be = 0.25;      % equivalent bandwidth
Nf = 1024;                          % 频点数
f = linspace(-0.5, 0.5, Nf);        % 归一化频率范围 [-0.5, 0.5]
df = f(2) - f(1);                  % 频率分辨率
Ex = 1e1;
Q_ub = 1;
c = 1500;

% target info
v_target = 5; % m/s
eta = 1 + 2*v_target / c;
fd = bandwidth * (eta - 1) / bandwidth;
fd_num = round(fd / df);

% CM sequence generation para
max_iter = 5000;
tolerance = 1e-1;
M = 1000;  % 时域信号的点数


% ------------------------------------------------------------------------%
[Hf, Pc, Pn, H2] = environment_init(Nf, bandwidth, f);
Hf_norm = Hf / max(Hf); H2_norm = abs(Hf_norm).^2;                 % 归一化
Pc_norm = Pc / max(Pc);                 % 归一化
Pn_norm = Pn / max(Pn);

% visualize |H(f)|^2, Pc, Pn
% figure;
% plot(f, H2, 'k--', 'LineWidth', 1.5); hold on;
% plot(f, Pc, 'g--', 'LineWidth', 1.5);
% plot(f, Pn, 'k--', 'LineWidth', 1);
% legend('Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)', 'Location', 'bestoutside');
% title('Target and Interference Spectrum');
% xlabel('Normalized Frequency');
% ylabel('ESD and PSD');
% grid on;
figure;
plot(f, H2_norm, 'k--', 'LineWidth', 1.5); hold on;
plot(f, Pc_norm, 'g--', 'LineWidth', 1.5);
plot(f, Pn_norm, 'k--', 'LineWidth', 1);
legend('Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)', 'Location', 'bestoutside');
title('Target and Interference Spectrum');
xlabel('Frequency');
ylabel('ESD and PSD');
grid on;

% methods_to_run = {'waterfill', 'joint', 'QfuncComb', 'lfm', 'bangbang'};     % waterfill, bangbang, joint, Qfunc, QfuncComb， lfm
methods_to_run = {'QfuncComb'};     % waterfill, bangbang, joint, Qfunc, QfuncComb， lfm

results = struct();
for i = 1:length(methods_to_run)
    method = methods_to_run{i};
    fprintf('**** %s 算法 ****\n', method);

    switch method
        case 'waterfill'
            X_ESD = method_waterfill(H2, Pc, Pn, f, Ex);
        case 'bangbang'
            X_ESD = method_bangbang(H2, Pc, f, Ex);
        case 'joint'
            X_ESD = method_joint(H2, Pc, f, Ex, Be);
        case 'Qfunc'
            X_ESD = method_Qfunc(H2, Pc, Pn, f, Ex, Be, Q_ub, fd_num);
        case 'QfuncComb'
            X_ESD = method_Qfunc_comb(H2, Pc, Pn, f, Ex, fd);
        case 'lfm'
            B_lfm = 1;  % LFM 带宽（可根据需求调整），归一化值 ∈ [0, 1]
            X_ESD = method_lfm(f, Ex, B_lfm);

        otherwise
            warning('未知方法: %s', method);
            continue;
    end
    results.(method).X_ESD = X_ESD;

    results.(method).SINR = SINR_compute(X_ESD, H2, Pc, Pn);
    fprintf('SINR(dB): %.4f\n', 10 * log10(results.(method).SINR));
    fprintf('fd处Q函数值为： %.4f \n',  Qfunc(X_ESD, fd_num, f, Ex));
    fprintf('注水算法等效带宽： %.4f \n', equivalent_bandwidth(X_ESD, f));
    fprintf('能量： %.4f \n', energy(X_ESD, f))

    [signal, ESD_synthesized] = synthesize_signal_from_ESD(X_ESD.', M, max_iter, tolerance);
    results.(method).signal = signal;
    results.(method).autocorr = xcorr(signal, 'normalized'); % 归一化自相关
    results.(method).lags = -length(signal)+1 : length(signal)-1; % 时延序列
    results.(method).ESD_synthesized = ESD_synthesized;

end

% 环境PSD+波形ESD的归一化频谱
figure;
plot(f, H2_norm, 'k--', 'LineWidth', 0.5); hold on;
plot(f, Pc_norm, 'g--', 'LineWidth', 0.5);
plot(f, Pn_norm, 'k--', 'LineWidth', 1);
colors = {'b','r', 'm', 'c', 'g', 'y'};  % 自动配色用
for i = 1:length(methods_to_run)
    method = methods_to_run{i};
    X_ESD_norm = normalize_esd(results.(method).X_ESD);
    plot(f, X_ESD_norm, '-', 'LineWidth', 0.5, 'Color', colors{i});
end
legend_labels = [{'Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)'}, ...
                 strcat(methods_to_run, ' ESD')];
legend(legend_labels, 'Location', 'bestoutside');
title('Target and Interference Spectrum');
xlabel('Normalized Frequency');
ylabel('Normalized ESD and PSD');
grid on;

% 各自方法波形的ESD
figure;
for i = 1:length(methods_to_run)
    method = methods_to_run{i};
    plot(f, results.(method).X_ESD, 'LineWidth', 0.5, 'Color', colors{i}); hold on;
end
legend(methods_to_run, 'Location', 'bestoutside');
title('Waveform Power Spectrum |X(f)|^2');
xlabel('Normalized Frequency');
ylabel('ESD |X(f)|^2');
grid on;

% 时域波形自相关（AF 0多普勒切面）
figure;
for i = 1:length(methods_to_run)
    method = methods_to_run{i};
    ac = results.(method).autocorr;
    lags = results.(method).lags;
    subplot(2, 3, i);
    plot(lags, 20*log10(abs(ac)), 'LineWidth', 0.5, 'Color', colors{i});
    title([strrep(method, '_', '\_') ' autocorr']);
    % legend(strcat(methods_to_run, ' autocorr'), 'Location', 'bestoutside');
    xlabel('时延（采样点）');
    ylabel('归一化自相关幅值 (dB)');
    title('合成信号的自相关函数（对数刻度）');
    ylim([-60, 0]);
    grid on;
end


%% -------------------------测试部分-----------------------------------------------%

% 比对给定的ESD和恒模生成的ESD
method = 'QfuncComb';
figure;
plot(f, normalize_esd(results.(method).X_ESD), 'b--', 'LineWidth', 0.5);
hold on;
plot(f, results.(method).ESD_synthesized, 'r', 'LineWidth', 0.5);
xlabel('频率');
ylabel('ESD');
title('目标ESD与合成ESD对比');
legend('目标ESD', '合成ESD');
grid on;

% 恒模生成波形的时域
figure;
plot(real(results.(method).signal), 'b', 'LineWidth', 0.5);
% hold on;
% plot(imag(a_final), 'r--', 'LineWidth', 1.5);
xlabel('时间点');
ylabel('幅值');
title('合成时域信号');
legend('实部');
grid on;

% figure;
% titles = {'Original Water-fill', 'bangbang ESD', 'Joint Constraints ESD', 'Qfunc ESD'};
% data = {autocorr_WF, autocorr_bangbang, autocorr_joint, autocorr_Qfunc};
% colors = {'b', 'k', 'r', 'm'};
% for i = 1:4
%     subplot(2,2,i);
%     plot(lags, 20*log10(abs(data{i})), colors{i}, 'LineWidth', 0.5);
%     xlabel('时延（采样点）');
%     ylabel('归一化自相关幅值 (dB)');
%     title(titles{i});
%     ylim([-60, 0]);
%     grid on;
% end
% sgtitle('合成信号的自相关函数（对数刻度）');

% %% 测试节
% SINR_WF = sum(H2.* X_ESD_WF) / sum(X_ESD_WF .* Pc + Pn);
% SINR_comb = sum(H2.* X_ESD_comb) / sum(X_ESD_comb .* Pc + Pn);
% fprintf('Original WF SINR (dB): %.4f\n', 10 * log10(SINR_WF));
% fprintf('Original WF SINR (dB): %.4f\n', 10 * log10(SINR_comb));
%
% figure;
% plot(f, H2_norm, 'k--', 'LineWidth', 1.5); hold on;
% plot(f, Pc_norm, 'g--', 'LineWidth', 1.5);
% plot(f, Pn_norm, 'k--', 'LineWidth', 1);
% plot(f, X_ESD_WF_norm, 'b-', 'LineWidth', 1);
% plot(f, X_ESD_comb_norm, 'm-', 'LineWidth', 0.5);
% legend('Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)', 'Original WF ESD', 'Comb ESD', 'Location', 'bestoutside');
% title('Target and Interference Spectrum');
% xlabel('Normalized Frequency');
% ylabel('Nomalized ESD and PSD');
% grid on;

%% ------------------------  封装函数    ------------------------------------%
function [Hf, Pc, Pn, H2] = environment_init(Nf, bandwidth, f)
% target frequency response H(f),  consist of several Gaussian spectra. ESD = abs(Hf).^2
Hf = zeros(1, Nf);
% Data from Table 1, | center frequency | weight | variance | of each target Gaussian spectrum
target_params = [
    -0.42, 0.3, 1.70;
    -0.38, 0.5, 1.70;
    -0.25, 0.8, 1.70;
    -0.20, 1.0, 8.00;
    -0.07, 1.0, 8.00;
    0.08, 1.0, 1.70;
    0.15, 0.2, 1.70;
    0.28, 0.25, 1.70;
    0.31, 0.8, 8.00;
    0.39, 0.9, 8.00
    ];
target_params(:, 3) = target_params(:, 3) * 1e2 / bandwidth;
for i = 1:size(target_params, 1)
    Hf = Hf + target_params(i, 2) / (sqrt(2*pi)*target_params(i, 3)) ...
        * exp(-(f - target_params(i, 1)).^2 / (2 * target_params(i, 3).^2));
end
H2 = abs(Hf).^2;

% Clutter PSD
Pc = zeros(1, Nf);
clutter_params = [                 % 杂波参数（频点、权重、方差）
    -0.25, 1.0, 1.3;
    0.20, 1.0, 1.3
    ];
clutter_params(:, 3) = clutter_params(:, 3) * 1e3 / bandwidth;
for i = 1:size(clutter_params, 1)
    Pc = Pc + clutter_params(i, 2) / (sqrt(2*pi)*clutter_params(i, 3)) ...
        * exp(-(f - clutter_params(i, 1)).^2 / (2 * clutter_params(i, 3).^2));
end

% 噪声PSD
Pn = ones(1, Nf);                   % 白噪声（单位功率）

% 按照论文中TCR = CNR = 0 dB, 进行能量归一化
energy_target = trapz(f, H2);    % 目标的能量
energy_clutter = trapz(f, Pc);           % 杂波的能量
energy_noise = trapz(f, Pn);             % 噪声的能量
scale_target = 1 / energy_target;   % 目标的缩放因子
scale_clutter = 1 / energy_clutter; % 杂波的缩放因子
scale_noise = 1 / energy_noise;     % 噪声的缩放因子
Hf = Hf * sqrt(scale_target);
Pc = Pc * scale_clutter;
Pn = Pn * scale_noise;
H2 = abs(Hf).^2;
end

function X_ESD = method_waterfill(H2, Pc, Pn, f, Ex)
% Water filling (max eqution 1)
% Refer to:
% Theory and application of SNR and mutual information matched illumination waveforms
    df = f(2) - f(1);
    Bf = sqrt(H2 .* Pn) ./ Pc;
    Df = sqrt(Pn ./ H2);
    % 注水算法：通过二分法确定常数 A = 1 / sqrt(lambda)
    A_low = min(Df); % A 的下界
    A_high = max(Df + Ex ./ Bf); % A 的上界（估计值）
    tolerance = 1e-6; % 容差
    max_iter = 1000; % 最大迭代次数

    for iter = 1:max_iter
        A = (A_low + A_high) / 2;
        X_ESD = max(0, Bf .* (A - Df));
        current_Ex = trapz(f, X_ESD); % trapz：积分，计算当前能量
        if abs(current_Ex - Ex) < tolerance
            fprintf('注水算法在 %d 次迭代后收敛。\n', iter);
            break; % 满足能量约束
        elseif current_Ex < Ex
            A_low = A; % 能量不足，增加 A
        else
            A_high = A; % 能量过多，减小 A
        end
    end
end

function X_ESD = method_bangbang(H2, Pc, f, Ex)
    df = f(2) - f(1);
    Nf = length(f);
    cvx_clear
    cvx_begin quiet
        variable X_ESD(1, Nf) nonnegative;
        minimize(sum((Pc - H2) .* X_ESD));   % 目标函数
        subject to
            sum(X_ESD) * df <= Ex;           % 发射能量约束
    cvx_end
end

function X_ESD = method_joint(H2, Pc, f, Ex, Be)
% corrected SINR (max eqution 2, Joint Constraints Waveform Spectrum Design)
    df = f(2) - f(1);
    Nf = length(f);
    cvx_clear
    cvx_begin quiet
        variable X_ESD(1, Nf) nonnegative;      % 功率谱密度 |X(f)|^2
        minimize(sum((Pc - H2) .* X_ESD));      % 目标函数
        subject to
            sum(X_ESD) * df == Ex;              % 能量约束
            sum(X_ESD.^2) * df <= Ex^2 / Be;    % 等效带宽约束 (Aτ ≤ 1/Be)
    cvx_end
end

function X_ESD = method_Qfunc(H2, Pc, Pn, f, Ex, Be, Q_ub, fd_num)
% corrected SINR (max eqution 2, Q function Constraints Waveform Spectrum Design)
    df = f(2) - f(1);
    Nf = length(f);
    Aeq = ones(1, Nf);
    beq = Ex / df;
    x0 = (Ex / (Nf * df)) * ones(Nf, 1);

    cvx_clear
    fun = @(x) (Pc - H2) * x;     % corrected SINR
    % fun = @(x) -sum((H2' .* x) ./ (Pc' .* x + Pn')) * df;       % SINR Upper bound
    % 能量约束 (线性等式约束)
    Aeq = ones(1, Nf);
    beq = Ex / df;
    % 非线性约束
    nonlcon = @(x) deal([
        x' * x * df - (2* Ex^2 / Be);                          % 等效带宽约束 (x^2项)
        x' * circshift(x, fd_num) * df - Q_ub * Ex^2        % 模糊函数约束
        x' * circshift(x, fd_num-1) * df - Q_ub * Ex^2        % 模糊函数约束
        % x' * circshift(x, fd_num+1) * df - Q_ub * Ex^2        % 模糊函数约束
        ], []);
    % 设置 fmincon 选项
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...     % 'interior-point' 'sqp'
        'Display', 'final', ...
        'MaxIterations', 2000, ...
        'MaxFunctionEvaluations', 1e6);
    problem = createOptimProblem('fmincon', ...
        'objective', fun, ...
        'x0', x0, ...
        'Aeq', Aeq, ...
        'beq', beq, ...
        'lb', zeros(Nf, 1), ...
        'nonlcon', nonlcon, ...
        'options', options);
    % 设置 MultiStart 搜索器，执行多起点搜索，尝试 20 个随机初值
    ms = MultiStart('UseParallel', false, 'Display', 'none');  % 可设置 true 并行加速
    [x_best, fval_best, exitflag, output, manymins] = run(ms, problem, 1);
    X_ESD = x_best';
end

function X_ESD = method_Qfunc_comb(H2, Pc, Pn, f, Ex, fd)
    % corrected SINR Q func + Comb waveform
    Nf = length(f);
    X_ESD_WF = method_waterfill(H2, Pc, Pn, f, Ex);

    GC_r = 1.1;
    delta = 3 * fd; 
    f_start = -0.5;      % 起始频率

    comb_freqs = f_start;   % 第一个频点
    f_next = f_start + delta;
    while f_next <= 0.5
        comb_freqs(end+1) = f_next; %#ok<AGROW>
        delta = delta * GC_r;          % 间隔按等比增长
        f_next = f_next + delta;
    end

    X_ESD_comb = zeros(1, Nf);
    for k = 1:length(comb_freqs)
        [~, idx] = min(abs(f - comb_freqs(k))); % 找到最近频点
        X_ESD_comb(idx) = X_ESD_WF(idx);        % 使用注水niy结果幅值
    end
    % 能量归一化 
    X_ESD = X_ESD_comb * (trapz(f, X_ESD_WF) / trapz(f, X_ESD_comb));

end

function X_ESD = method_lfm(f, Ex, B_lfm)
    % f: 归一化频率向量 [-0.5, 0.5]
    % Ex: 总发射能量
    % B_lfm: LFM 波形的归一化带宽，例如0.2表示占频带的20%

    df = f(2) - f(1);
    Nf = length(f);
    X_ESD = zeros(1, Nf);

    f_low = -B_lfm / 2;
    f_high = B_lfm / 2;
    idx_range = (f >= f_low) & (f <= f_high);
    X_ESD(idx_range) = 1;

    % 能量归一化
    current_energy = sum(X_ESD) * df;
    X_ESD = X_ESD * (Ex / current_energy);
end

function Energy = energy(ESD, f)
    df = f(2) - f(1);
    Energy = sum(ESD) * df;
end

function equal_bw = equivalent_bandwidth(ESD, f)
    df = f(2) - f(1);
    equal_bw = (sum(ESD) * df)^2 / (sum(ESD.^2) * df);
end

function Q_val = Qfunc(ESD, fd_num, f, Ex)
    df = f(2) - f(1);
    Q_val = sum(ESD .* circshift(ESD, [0, fd_num])) * df;
end

function sinr_val = SINR_compute(X_ESD, H2, Pc, Pn)
    sinr_val = sum(H2 .* X_ESD) / sum(X_ESD .* Pc + Pn);
end

function ESD_norm = normalize_esd(ESD)
    ESD_norm = ESD / max(ESD);
end



