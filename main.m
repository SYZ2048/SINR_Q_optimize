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

bandwidth = 50e6; % 50 MHz
pulse_width = 2.5e-6; % 2.5 µs
Be = 0.25;
Nf = 1024;                          % 频点数
f = linspace(-0.5, 0.5, Nf);        % 归一化频率范围 [-0.5, 0.5]
df = f(2) - f(1);                  % 频率分辨率
Ex = 1e-2;

% ------------------------------------------------------------------------%
% target frequency response H(f),  consist of several Gaussian spectra. ESD = abs(Hf).^2
Hf = zeros(1, Nf);
% Data from Table 1, | center frequency | weight | variance | of each target Gaussian spectrum
target_params = [                   
    -0.42, 0.3, 1.70e6;
    -0.38, 0.5, 1.70e6;
    -0.25, 0.8, 1.70e6;
    -0.20, 1.0, 8.00e5;
    -0.07, 1.0, 8.00e5;
     0.08, 1.0, 1.70e6;
     0.15, 0.2, 1.70e6;
     0.28, 0.25, 1.70e6;
     0.31, 0.8, 8.00e5;
     0.39, 0.9, 8.00e5
];
target_params(:, 3) = target_params(:, 3) / bandwidth;
for i = 1:size(target_params, 1)
    Hf = Hf + target_params(i, 2) / (sqrt(2*pi)*target_params(i, 3)) ...
        * exp(-(f - target_params(i, 1)).^2 / (2 * target_params(i, 3).^2));
end
H2 = abs(Hf).^2;


% Clutter PSD
Pc = zeros(1, Nf);
clutter_params = [                 % 杂波参数（频点、权重、方差）
    -0.25, 1.0, 8.3e6;
     0.20, 1.0, 8.3e6
];
clutter_params(:, 3) = clutter_params(:, 3) / bandwidth;
for i = 1:size(clutter_params, 1)
    Pc = Pc + clutter_params(i, 2) / (sqrt(2*pi)*clutter_params(i, 3)) ...
        * exp(-(f - clutter_params(i, 1)).^2 / (2 * clutter_params(i, 3).^2));
end

% Pc = abs(Pc).^2;      %%%% not sure!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 噪声PSD
Pn = ones(1, Nf);                   % 白噪声（单位功率）

% 按照论文中TCR = CNR = 0 dB, 进行能量归一化
energy_target = trapz(f, H2);    % 目标的能量
energy_clutter = trapz(f, Pc);           % 杂波的能量
energy_noise = trapz(f, Pn);             % 噪声的能量
scale_target = 1 / energy_target;   % 目标的缩放因子
scale_clutter = 1 / energy_clutter; % 杂波的缩放因子
scale_noise = 1 / energy_noise;     % 噪声的缩放因子

% 对目标、杂波和噪声进行缩放，使它们的能量相等
Hf = Hf * sqrt(scale_target);
Pc = Pc * scale_clutter;
Pn = Pn * scale_noise;
H2 = abs(Hf).^2;

Hf_norm = Hf / max(Hf); H2_norm = abs(Hf_norm).^2;                 % 归一化
Pc_norm = Pc / max(Pc);                 % 归一化
Pn_norm = Pn / max(Pn);

%% ------------------------------------------------------------------------% 
% Water filling (max eqution 1)
% Refer to:
% Theory and application of SNR and mutual information matched illumination waveforms
%
Bf = sqrt(H2 .* Pn) ./ Pc;
Df = sqrt(Pn ./ H2);
% 注水算法：通过二分法确定常数 A = 1 / sqrt(lambda)
A_low = min(Df); % A 的下界
A_high = max(Df + Ex ./ Bf); % A 的上界（估计值）
tolerance = 1e-6; % 容差
max_iter = 1000; % 最大迭代次数

for iter = 1:max_iter
    A = (A_low + A_high) / 2;
    X_ESD_WF = max(0, Bf .* (A - Df));
    current_Ex = trapz(f, X_ESD_WF); % trapz：积分，计算当前能量

    if abs(current_Ex - Ex) < tolerance
        disp(iter);
        break; % 满足能量约束
    elseif current_Ex < Ex
        A_low = A; % 能量不足，增加 A
    else
        A_high = A; % 能量过多，减小 A
    end
end

X_ESD_WF_norm = X_ESD_WF / max(X_ESD_WF);

(trapz(f, X_ESD_WF))^2 / trapz(f, X_ESD_WF .^ 2)    % e equivalent bandwidth

% ------------------------------------------------------------------------% 
% corrected SINR (max eqution 2)
cvx_clear
cvx_begin
    variable X_ESD_bangbang(1,Nf) nonnegative;    % 功率谱密度 |X(f)|^2
    minimize(sum((Pc - H2) .* X_ESD_bangbang)); % 目标函数
    subject to
        sum(X_ESD_bangbang) * df <= Ex;         % 能量约束
        % sum(X_ESD_bangbang.^2) * df <= (Ex^2 / Be); % 等效带宽约束 (Aτ ≤ 1/Be)
cvx_end
X_ESD_bangbang_norm = X_ESD_bangbang / max(X_ESD_bangbang);

(trapz(f, X_ESD_bangbang))^2 / trapz(f, X_ESD_bangbang .^ 2)

% ------------------------------------------------------------------------% 
% corrected SINR (max eqution 2, Joint Constraints Waveform Spectrum Design)
cvx_clear
cvx_begin
    variable X_ESD_joint(1,Nf) nonnegative;    % 功率谱密度 |X(f)|^2
    minimize(sum((Pc - H2) .* X_ESD_joint)); % 目标函数
    subject to
        sum(X_ESD_joint) * df <= Ex;         % 能量约束
        sum(X_ESD_joint.^2) * df <= (Ex^2 / Be); % 等效带宽约束 (Aτ ≤ 1/Be)
cvx_end
X_ESD_joint_norm = X_ESD_joint / max(X_ESD_joint);

(trapz(f, X_ESD_joint))^2 / trapz(f, X_ESD_joint .^ 2)

% According to Algorithm 1, NOT FINISHED!!!!
% X_ESD_Joint = ones(1, Nf) * Ex / (Nf * df); % 初始均匀分布
% v = 0;                            % Lagrange乘数初始化
% lambda = 0;                       % Lagrange乘数初始化
% step = 0.01;
% % 迭代求解KKT条件
% for iter = 1:max_iter
%     x_prev = X_ESD_Joint;
%     for i = 1:Nf
%         numerator = -Pc(i) + ESD(i) - v;
%         if numerator > 0
%             X_ESD_Joint(i) = numerator / (2 * lambda);
%         else
%             X_ESD_Joint(i) = 0;
%         end
%     end
% 
%     % 更新v和lambda（Algorithm 1）
%     v = v + step; % 搜索步长
%     lambda = (sum(X_ESD_Joint) - Ex/df) / (2 * Ex/df);
% 
%     % 检查收敛条件
%     if norm(X_ESD_Joint - x_prev) < tolerance
%         break;
%     end
% end
% X_ESD_Joint = X_ESD_Joint / sum(X_ESD_Joint) * (Ex / df);
% ------------------------------------------------------------------------% 
% corrected SINR (max eqution 2, Q function constrain)
% cvx_clear
% cvx_begin
%     variable X_ESD_joint(1,Nf) nonnegative;    % 功率谱密度 |X(f)|^2
%     minimize(sum((Pc - H2) .* X_ESD_joint)); % 目标函数
%     subject to
%         sum(X_ESD_joint) * df <= Ex;         % 能量约束
%         sum(X_ESD_joint.^2) * df <= (Ex^2 / Be); % Q function constrain
% cvx_end
% X_ESD_joint_norm = X_ESD_joint / max(X_ESD_joint);




%% 绘制结果
figure;
plot(f, H2_norm, 'k--', 'LineWidth', 1.5); hold on;
plot(f, Pc_norm, 'g--', 'LineWidth', 1.5);
plot(f, Pn_norm, 'k--', 'LineWidth', 1);
plot(f, X_ESD_WF_norm, 'b-', 'LineWidth', 1);
plot(f, X_ESD_bangbang_norm, 'r-', 'LineWidth', 1.5);
plot(f, X_ESD_joint_norm, 'y-', 'LineWidth', 1);
legend('Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)', 'Original WF ESD', 'bangbang ESD', 'joint ESD', 'Location', 'bestoutside');
title('Target and Interference Spectrum');
xlabel('Normalized Frequency');
ylabel('Nomalized ESD and PSD');
grid on;

% 不归一化啥也看不见了……
figure;
plot(f, H2, 'k--', 'LineWidth', 1.5); hold on;
plot(f, Pc, 'g--', 'LineWidth', 1.5);
plot(f, Pn, 'k--', 'LineWidth', 1);
plot(f, X_ESD_WF, 'b-', 'LineWidth', 1);
plot(f, X_ESD_bangbang, 'r-', 'LineWidth', 1.5);
plot(f, X_ESD_joint, 'y-', 'LineWidth', 1);
legend('Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)', 'Original WF ESD', 'bangbang ESD', 'joint ESD', 'Location', 'bestoutside');
title('Target and Interference Spectrum');
xlabel('Normalized Frequency');
ylabel('ESD and PSD');
grid on;


% figure;
% plot(f, X_wf, 'g', 'LineWidth', 1.5); hold on;
% plot(f, X_corrected, 'm', 'LineWidth', 1.5);
% legend('Original WF', 'Corrected Metric');
% title('Waveform Spectrum Design');
% xlabel('Normalized Frequency');
% ylabel('Energy Spectral Density');


%% 计算SINR
SINR_WF = sum(abs(Hf).^2 .* X_ESD_WF) / sum(X_ESD_WF .* Pc + Pn);
SINR_bangbang = sum(abs(Hf).^2 .* X_ESD_bangbang) / sum(X_ESD_bangbang .* Pc + Pn);
X_ESD_joint = sum(abs(Hf).^2 .* X_ESD_joint) / sum(X_ESD_joint .* Pc + Pn);

% SINR_corrected = sum(abs(Hf).^2 .* X_corrected) / sum(X_corrected .* Pc + Pn);
fprintf('Original WF SINR: %.4f\n', SINR_WF);
fprintf('Corrected Metric SINR: %.4f\n', SINR_bangbang);
fprintf('Corrected Metric SINR: %.4f\n', X_ESD_joint);

%% 时域恒模序列合成（Algorithm 2）
max_iter = 2000;
tolerance = 1e-1;
M = 1000;  % 时域信号的点数
disp(1);
[signal_WF, ESD_WF_synthesized] = synthesize_signal_from_ESD(X_ESD_WF_norm.', M, max_iter, tolerance);
disp(2);
[signal_bangbang, ESD_bangbang_synthesized] = synthesize_signal_from_ESD(X_ESD_bangbang_norm.', M, max_iter, tolerance);
disp(3);
[signal_joint, ESD_joint_synthesized] = synthesize_signal_from_ESD(X_ESD_joint_norm.', M, max_iter, tolerance);
disp(4);
autocorr_WF = xcorr(signal_WF, 'normalized'); % 归一化自相关
autocorr_bangbang = xcorr(signal_bangbang, 'normalized'); % 归一化自相关
autocorr_joint = xcorr(signal_joint, 'normalized'); % 归一化自相关
lags = -length(signal_joint)+1 : length(signal_joint)-1; % 时延序列

figure;
plot(lags, 20*log10(abs(autocorr_WF)), 'b', 'LineWidth', 0.5);hold on;
plot(lags, 20*log10(abs(autocorr_bangbang)), 'k', 'LineWidth', 0.5);
plot(lags, 20*log10(abs(autocorr_joint)), 'r', 'LineWidth', 0.5);
xlabel('时延（采样点）');
ylabel('归一化自相关幅值 (dB)');
title('合成信号的自相关函数（对数刻度）');
ylim([-60, 0]); % 限制动态范围以观察旁瓣
grid on;






