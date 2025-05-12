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
Be = 0.1;
Nf = 1024;                          % 频点数
f = linspace(-0.5, 0.5, Nf);        % 归一化频率范围 [-0.5, 0.5]
df = f(2) - f(1);                  % 频率分辨率
Ex = 1e-2;

% ------------------------------------------------------------------------%
% target frequency response H(f),  consist of several Gaussian spectra. ESD = abs(Hf).^2
Hf = zeros(1, Nf);
% Data from Table 1, center frequency, weight and variance of each target Gaussian spectrum
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
Hf_norm = Hf / max(Hf);                 % 归一化
H2 = abs(Hf).^2;
H2_norm = abs(Hf_norm).^2;

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
%%%% not sure!!!!!!!!!!
Pc = abs(Pc).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pc_norm = Pc / max(Pc);                 % 归一化

% 噪声PSD
Pn = ones(1, Nf);                   % 白噪声（单位功率）

%% ------------------------------------------------------------------------% 
% Water filling (max eqution 1)
% Refer to:
% Theory and application of SNR and mutual information matched illumination waveforms
%
Bf = sqrt(abs(Hf).^2 .* Pn) ./ Pc;
Df = sqrt(Pn ./ abs(Hf).^2);
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
        break; % 满足能量约束
    elseif current_Ex < Ex
        A_low = A; % 能量不足，增加 A
    else
        A_high = A; % 能量过多，减小 A
    end
end

X_ESD_WF_norm = X_ESD_WF / max(X_ESD_WF);
X_spectrum_WF = sqrt(X_ESD_WF) .* exp(1i * angle(Hf)); % 保持相位信息

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
plot(f, Pn, 'k--', 'LineWidth', 1);
plot(f, X_ESD_WF_norm, 'b-', 'LineWidth', 1);
plot(f, X_ESD_bangbang_norm, 'r-', 'LineWidth', 1.5);
plot(f, X_ESD_joint_norm, 'y-', 'LineWidth', 1);
legend('Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)', 'Original WF ESD', 'bangbang ESD', 'joint ESD', 'Location', 'bestoutside');
title('Target and Interference Spectrum');
xlabel('Normalized Frequency');
ylabel('Nomalized ESD and PSD');
grid on;

% 不归一化啥也看不见了……
% figure;
% plot(f, H2, 'k--', 'LineWidth', 1.5); hold on;
% plot(f, Pc, 'g--', 'LineWidth', 1.5);
% plot(f, Pn, 'k--', 'LineWidth', 1);
% plot(f, X_ESD_WF, 'b-', 'LineWidth', 1);
% plot(f, X_ESD_bangbang, 'r-', 'LineWidth', 1.5);
% plot(f, X_ESD_joint, 'y-', 'LineWidth', 1);
% legend('Target ESD |H(f)|^2', 'Clutter PSD P_c(f)', 'Noise PSD P_n(f)', 'Original WF ESD', 'bangbang ESD', 'joint ESD', 'Location', 'bestoutside');
% title('Target and Interference Spectrum');
% xlabel('Normalized Frequency');
% ylabel('ESD and PSD');
% grid on;


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
% X_spectrum = sqrt(x) .* exp(1i * 2*pi*rand(1,Nf)); % 随机相位
% y = ifft(ifftshift(X_spectrum));    % IFFT得到时域信号
% y = y / max(abs(y));               % 归一化幅值
% 
% figure;
% plot(abs(y), 'LineWidth', 1.5);
% title('Time-Domain CM Sequence');
% xlabel('Time Samples');
% ylabel('Amplitude');








