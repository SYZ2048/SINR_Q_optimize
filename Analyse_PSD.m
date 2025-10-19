% ------------------------------------------------------------------------%
% Analyse PSD and |H^2(f)| from Real Data
% SYZ-20250707
% ------------------------------------------------------------------------%
%% 加载数据
clc;clear;
close all;

% noise_only - 纯背景噪声数据
% clutter_noise - 背景噪声+混响数据(无目标)
% target_clutter_noise - 背景噪声+混响+目标回波数据(有目标)
% -------------------------------------------------------------------------
%                   无目标     有目标     无目标放大25倍    有目标放大25倍
%   20mvrms           7         13            6               12
%   100               3         10            2               11
%   500               4         9             5                8
% -------------------------------------------------------------------------
test_num = 2;
if test_num == 1
    noise_only = tdmsread("./20251014/LFM100k_B40k_fs400k_200ms_20ms_1.tdms");
    clutter_noise = tdmsread("./20251014/LFM100k_B40k_fs400k_200ms_20ms_4.tdms");
    target_noise = tdmsread("./20251014/LFM100k_B40k_fs400k_200ms_20ms_9.tdms");
    T_cycle = 200e-3;
    Tp = 20e-3;
    fc = 100e3;
    Bw = 40e3;
elseif test_num == 2
    noise_only = tdmsread("./20251015/LFM60k_B40k_fs400k_100ms_5ms_1.tdms");
    clutter_noise = tdmsread("./20251015/LFM60k_B40k_fs400k_100ms_5ms_14.tdms");
    target_noise = tdmsread("./20251015/LFM60k_B40k_fs400k_100ms_5ms_11.tdms");
    T_cycle = 100e-3;
    Tp = 5e-3;
    fc = 60e3;
    Bw = 40e3;
elseif test_num == 3
    noise_only = tdmsread("./20251015/LFM60k_B40k_fs400k_100ms_5ms_1.tdms");
    clutter_noise = tdmsread("./20251015/LFM80k_B40k_fs400k_100ms_5ms_13.tdms");
    target_noise = tdmsread("./20251015/LFM80k_B40k_fs400k_100ms_5ms_12.tdms");
    T_cycle = 100e-3;
    Tp = 5e-3;
    fc = 80e3;
    Bw = 40e3;
elseif test_num == 4
    noise_only = tdmsread("./20251016/noise.tdms");
    clutter_noise = tdmsread("./20251016/LFM10k-20k_fs400k_100ms_5ms_5.tdms");
    target_noise = tdmsread("./20251016/LFM10k-20k_fs400k_100ms_5ms_4.tdms");
    T_cycle = 100e-3;
    Tp = 5e-3;
    fc = 12.5e3;
    Bw = 15e3;
elseif test_num == 5
    noise_only = tdmsread("./20251016/noise.tdms");
    clutter_noise = tdmsread("./20251016/LFM10k-20k_fs400k_500ms_200ms_2.tdms");
    target_noise = tdmsread("./20251016/LFM10k-20k_fs400k_500ms_200ms_3.tdms");
    T_cycle = 500e-3;
    Tp = 200e-3;
    fc = 12.5e3;
    Bw = 15e3;
end


noise_only_chanvals = noise_only{1:1}{:,:};  % N * 2
noise_only_recv_sig = noise_only_chanvals(:, 1);
noise_only_window = noise_only_chanvals(:, 2);

clutter_noise_chanvals = clutter_noise{1:1}{:,:};  % N * 2
clutter_noise_recv_sig = clutter_noise_chanvals(:, 1);
clutter_noise_window = clutter_noise_chanvals(:, 2);

target_noise_chanvals = target_noise{1:1}{:,:};  % N * 2
target_noise_recv_sig = target_noise_chanvals(:, 1);
target_noise_window = target_noise_chanvals(:, 2);

% 参数设置
fs = 400e3; % 采样率
fmin = fc - Bw/2;
fmax = fc + Bw/2;
Np = fs * Tp;
N_cycle = fs * T_cycle;
c_v = 1500;

analyse_cycle = 10;
N_analyse = N_cycle * analyse_cycle;
t_mf = (0 : N_analyse-1) / fs;

threshold = 2;  % Reference signal values vary in [0, 3]
[start_idx_noise, end_idx_noise] = ...
    findAnalysisInterval(noise_only_chanvals, threshold, N_analyse);
noise_only_recv_sig = noise_only_recv_sig(start_idx_noise:end_idx_noise);
noise_only_window = noise_only_window(start_idx_noise:end_idx_noise);
% figure();
% subplot(2,1,1);plot(t_mf*1e3, noise_only_recv_sig);title('Noise - Receive Signal in Time Domain - Selected Cycles');xlabel('Time(ms)')
% subplot(2,1,2);plot(t_mf*1e3, noise_only_window);title('Noise - Receive Signal Window in Time Domain - Selected Cycles');xlabel('Time(ms)')

[start_idx_clutter, end_idx_clutter] = ...
    findAnalysisInterval(clutter_noise_chanvals, threshold, N_analyse);
clutter_noise_recv_sig = clutter_noise_recv_sig(start_idx_clutter:end_idx_clutter);
clutter_noise_window = clutter_noise_window(start_idx_clutter:end_idx_clutter);
% figure();
% subplot(2,1,1);plot(t_mf*1e3, clutter_noise_recv_sig);title('Clutter - Receive Signal in Time Domain - Selected Cycles');xlabel('Time(ms)')
% subplot(2,1,2);plot(t_mf*1e3, clutter_noise_window);title('Clutter - Receive Signal Window in Time Domain - Selected Cycles');xlabel('Time(ms)')

[start_idx_target, end_idx_target] = ...
    findAnalysisInterval(target_noise_chanvals, threshold, N_analyse);
target_noise_recv_sig = target_noise_recv_sig(start_idx_target:end_idx_target);
target_noise_window = target_noise_window(start_idx_target:end_idx_target);
% figure();
% subplot(2,1,1);plot(t_mf*1e3, target_noise_recv_sig);title('Target - Receive Signal in Time Domain - Selected Cycles');xlabel('Time(ms)')
% subplot(2,1,2);plot(t_mf*1e3, target_noise_window);title('Target - Receive Signal Window in Time Domain - Selected Cycles');xlabel('Time(ms)')

%% Matched Filter: 
k = Bw/Tp;
t = (0 : Np-1) / fs;
lfm = sin(2*pi*fmin*t +pi*k*t.^2);
% figure();plot(lfm);title('Transmit Signal');
% spectrogram(lfm,256,128,100,fs,'yaxis');
% transmit_signal = target_noise_recv_sig(start_idx_target + 300:start_idx_target + 8300);
transmit_signal = lfm;

smf_target = matchedfilter(target_noise_recv_sig,transmit_signal,fs);
smf_target = smf_target / max(abs(smf_target));
% figure();plot(t_mf*1e3, smf_target);title('Target - Matched Filter Result- Selected Cycles');xlabel('Time(ms)')
% xlim([100 120]);ylim([-0.01 0.01]);

smf_clutter = matchedfilter(clutter_noise_recv_sig,transmit_signal,fs);
smf_clutter = smf_clutter / max(abs(smf_clutter));
% figure();plot(t_mf*1e3, smf_clutter);title('Clutter - Matched Filter Result- Selected Cycles');xlabel('Time(ms)')
% xlim([100 120]);ylim([-0.01 0.01]);

%% 计算功率谱密度和目标频率响应
delay = 10.285e-3;
delay_direct_blast = 1.06e-3;
delay_inference = 11.4e-3;
distance = delay * c_v / 2

nfft = 512; % FFT点数
window = hann(nfft); % 窗函数
noverlap = nfft/2; % 重叠点数

% 计算功率谱密度(PSD)，结果为单边谱，P和f维度为(nfft/2+1) * 1 
length(noise_only_recv_sig(delay*fs:delay*fs + Np))

P_noise_avg = zeros(nfft/2+1, 1);   % N * 1
P_clutter_noise_avg = zeros(nfft/2+1, 1);
P_target_clutter_noise_avg = zeros(nfft/2+1, 1);

% --pwelch进行功率谱密度估计
%       对于实数信号输入，默认返回单边onesided功率谱，范围[0, fs/2]
%       对于复数信号输入，默认值为 'twosided'，范围为[0,fs]
%       'centered' - 对于实数值或复数值输入 x，范围为 [–fs/2, fs/2] 
[P_tx, f] = pwelch(transmit_signal, window, noverlap, nfft, fs);
num_cycles = analyse_cycle-1;
for cycle = 0:num_cycles-1
    cycle_start = cycle * N_cycle + delay*fs;
    
    [P_noise, f] = pwelch(noise_only_recv_sig(cycle_start:cycle_start + Np), window, noverlap, nfft, fs);
    [P_clutter_noise, f] = pwelch(clutter_noise_recv_sig(cycle_start:cycle_start + Np), window, noverlap, nfft, fs);
    [P_target_clutter_noise, f] = pwelch(target_noise_recv_sig(cycle_start:cycle_start + Np), window, noverlap, nfft, fs);
    

    P_noise_avg = P_noise_avg + P_noise;
    P_clutter_noise_avg = P_clutter_noise_avg + P_clutter_noise;
    P_target_clutter_noise_avg = P_target_clutter_noise_avg + P_target_clutter_noise;

end
P_noise_avg = P_noise_avg / num_cycles;
P_clutter_noise_avg = P_clutter_noise_avg / num_cycles;
P_target_clutter_noise_avg = P_target_clutter_noise_avg / num_cycles;
df = f(2)-f(1);
% -------------- Test Start ----------------------
% 多周期平均PSD可视化
figure;
plot(f/1e3, 10*log10(P_noise_avg), 'b', 'LineWidth', 1.2); hold on;
plot(f/1e3, 10*log10(P_clutter_noise_avg), 'r', 'LineWidth', 1.2);hold on;
plot(f/1e3, 10*log10(P_target_clutter_noise_avg), 'g', 'LineWidth', 1.2); hold on;
% plot(f/1e3, 10*log10(P_tx), 'm', 'LineWidth', 1.2);
title('多周期平均PSD');
xlabel('频率 (kHz)');
ylabel('PSD dB');
legend('Noise', 'Clutter+Noise', 'Clutter+Noise+|H(f)|^2', 'Location', 'best');
grid on;

[P_direct_blast, f] = pwelch(target_noise_recv_sig(floor(delay_direct_blast*fs):floor(delay_direct_blast*fs) + Np), window, noverlap, nfft, fs);
figure;
plot(f/1e3, 10*log10(P_direct_blast));
title('直达波PSD');
xlabel('频率 (kHz)');
ylabel('PSD dB');
grid on;

% 查看单周期的噪声PSD，混响PSD，目标|H(f)^2|，发射波形|S(f)^2|
% figure;
% plot(f/1e3, 10*log10(P_noise), 'b', 'LineWidth', 1.2); hold on;
% plot(f/1e3, 10*log10(P_clutter_noise), 'r', 'LineWidth', 1.2);hold on;
% plot(f/1e3, 10*log10(P_target_clutter_noise), 'g', 'LineWidth', 1.2); hold on;
% % plot(f/1e3, 10*log10(P_tx), 'm', 'LineWidth', 1.2);
% % title('发射波形 S|f|^2');
% title('')
% xlabel('频率 (kHz)');
% ylabel('PSD dB');
% legend('Noise', 'Clutter+Noise', 'Clutter+Noise+|H(f)|^2', '|S(f)|^2',  'Location', 'best');
% grid on;

% 多周期FFT
P_fft_avg_target = zeros(nfft/2+1, 1);
P_fft_avg_direct = zeros(nfft/2+1, 1);
P_fft_avg_clutter = zeros(nfft/2+1, 1);
P_fft_avg_inference = zeros(nfft/2+1, 1);

for cycle = 0:num_cycles-1
    cycle_start = cycle * N_cycle;

    % 目标回波段
    target_echo = smf_target(cycle_start + delay*fs : cycle_start + delay*fs + Np);
    X_fft = fft(target_echo, nfft);
    P_fft = abs(X_fft).^2 / length(target_echo);
    P_fft_single = P_fft(1:nfft/2+1);
    P_fft_single(2:end-1) = 2 * P_fft_single(2:end-1);
    P_fft_avg_target = P_fft_avg_target + P_fft_single(:);  

    % 11.4ms干扰波段
    inference = smf_target(cycle_start + floor(delay_inference*fs) : cycle_start + floor(delay_inference*fs) + Np);
    inference_fft = fft(inference, nfft);
    Pinference_fft = abs(inference_fft).^2 / length(inference);
    Pinference_fft_single = Pinference_fft(1:nfft/2+1);
    Pinference_fft_single(2:end-1) = 2 * Pinference_fft_single(2:end-1);
    P_fft_avg_inference = P_fft_avg_inference + Pinference_fft_single(:); 

    % 直达波段
    direct_blast = smf_target(cycle_start + floor(delay_direct_blast*fs) : cycle_start + floor(delay_direct_blast*fs) + Np);
    D_fft = fft(direct_blast, nfft);
    Pd_fft = abs(D_fft).^2 / length(direct_blast);
    Pd_fft_single = Pd_fft(1:nfft/2+1);
    Pd_fft_single(2:end-1) = 2 * Pd_fft_single(2:end-1);
    P_fft_avg_direct = P_fft_avg_direct + Pd_fft_single(:); 

    % 干扰（混响噪声）波段
    clutter = smf_clutter(cycle_start + floor(delay*fs) : cycle_start + floor(delay*fs) + Np);
    clutter_fft = fft(clutter, nfft);
    Pclutter_fft = abs(clutter_fft).^2 / length(clutter);
    Pclutter_fft_single = Pclutter_fft(1:nfft/2+1);
    Pclutter_fft_single(2:end-1) = 2 * Pclutter_fft_single(2:end-1);
    P_fft_avg_clutter = P_fft_avg_clutter + Pclutter_fft_single(:); 
end
P_fft_avg_target = P_fft_avg_target / num_cycles;
P_fft_avg_direct = P_fft_avg_direct / num_cycles;
P_fft_avg_clutter = P_fft_avg_clutter / num_cycles;
P_fft_avg_inference = P_fft_avg_inference / num_cycles;

figure;
subplot(2,2,1);
plot(f/1e3, 10*log10(P_fft_avg_target), 'b');
title('目标回波 - 多周期平均 FFT');
xlabel('频率 (kHz)'); ylabel('幅度 (dB)');
xlim([20 140]);
grid on;

subplot(2,2,2);
plot(f/1e3, 10*log10(P_fft_avg_inference), 'r');
title('11.4ms干扰 - 多周期平均 FFT');
xlabel('频率 (kHz)'); ylabel('幅度 (dB)');
xlim([20 140]);
grid on;

subplot(2,2,3);
plot(f/1e3, 10*log10(P_fft_avg_clutter), 'r');
title('混响噪声 - 多周期平均 FFT');
xlabel('频率 (kHz)'); ylabel('幅度 (dB)');
xlim([20 140]);
grid on;

subplot(2,2,4);
plot(f/1e3, 10*log10(P_fft_avg_direct), 'r');
title('直达波 - 多周期平均 FFT');
xlabel('频率 (kHz)'); ylabel('幅度 (dB)');
xlim([20 140]);
grid on;

% 单周期fft
% direct_blast = smf(1:1 + Np);
% X_fft = fft(target_echo, nfft);
% P_fft = abs(X_fft).^2 / length(target_echo);
% P_fft_single = P_fft(1:nfft/2+1);  % 单边谱
% P_fft_single(2:end-1) = 2 * P_fft_single(2:end-1);
% D_fft  = fft(direct_blast, nfft);
% Pd_fft = abs(D_fft).^2 / length(direct_blast);
% Pd_fft_single = Pd_fft(1:nfft/2+1);
% Pd_fft_single(2:end-1) = 2 * Pd_fft_single(2:end-1);
% figure;
% subplot(2,1,1);
% plot(f(1:nfft/2+1)/1e3, 10*log10(P_fft_single), 'b');
% title('Target Echo - FFT 周期图');
% xlabel('频率 (kHz)'); ylabel('幅度 (dB)');
% grid on;
% subplot(2,1,2);
% plot(f(1:nfft/2+1)/1e3, 10*log10(Pd_fft_single), 'r');
% title('Direct Blast - FFT 周期图');
% xlabel('频率 (kHz)'); ylabel('幅度 (dB)');
% grid on;
% -------------- Test End ----------------------


%% 估计目标响应 |H(f)|^2 从目标+混响+噪声中减去混响+噪声
P_clutter = (P_clutter_noise - P_noise) ./ P_tx;
P_clutter(P_clutter <= 0) = eps;
P_target = (P_target_clutter_noise - P_clutter_noise) ./ P_tx;
P_target(P_target <= 0) = eps;

P_noise_dB = 10*log10(P_noise);
P_clutter_dB = 10*log10(P_clutter);
P_clutter_noise_dB = 10*log10(P_clutter_noise);
P_target_dB = 10*log10(P_target);
% 绘制结果
figure;
% 背景噪声PSD
subplot(3,1,1);
plot(f/1e3, P_noise_dB);
title('背景噪声PSD');
xlabel('频率 (kHz)');
ylabel('PSD dB');
grid on;
% 混响PSD
subplot(3,1,2);
plot(f/1e3, P_clutter_dB);
title('混响 PSD');
xlabel('频率 (kHz)');
ylabel('PSD dB');
grid on;
% 目标响应 |H(f)|^2
subplot(3,1,3);
plot(f/1e3, P_target_dB);
title('目标响应 |H(f)|^2');
xlabel('频率 (kHz)');
ylabel('|H(f)|^2 dB');
grid on;


% -------------- Test Start ----------------------
% 
% figure;
% plot(f/1e3, P_target, 'k-.', 'LineWidth', 1.5); hold on;
% plot(f/1e3, P_clutter, 'g--', 'LineWidth', 1.5); hold on;
% plot(f/1e3, P_noise, 'b', 'LineWidth', 1.5); 
% lgd = legend('目标响应|H(f)|^2', '混响PSD', '噪声PSD');
% lgd.FontSize = 14;                          % 图例字号
% lgd.Box      = 'off';                       % 去掉边框（可选）
% xlabel('频率 kHz', 'FontSize', 18);
% ylabel('幅度', 'FontSize', 18);
% set(gca, 'FontSize', 17);       % 坐标刻度字号
% grid on;
% -------------- Test End ----------------------
%% 波形设计部分
f_lowbound = fc - Bw/4;
f_highbound = fc + Bw/4;
f_idx_start = find(f >= f_lowbound, 1, 'first');
f_idx_end = find(f <= f_highbound, 1, 'last');
% Nf = f_idx_end - f_idx_start + 1;
% f = f(f_idx_start:f_idx_end);
% Pc = P_clutter(f_idx_start:f_idx_end).';
% Pn = P_noise(f_idx_start:f_idx_end).';
% H2 = P_target(f_idx_start:f_idx_end).';
P_clutter([1:numel(P_clutter)] < f_idx_start | [1:numel(P_clutter)] > f_idx_end) = 0;
P_noise([1:numel(P_noise)] < f_idx_start | [1:numel(P_noise)] > f_idx_end) = 0;
P_target([1:numel(P_target)] < f_idx_start | [1:numel(P_target)] > f_idx_end) = 0;

figure;
plot(f/1e3, P_target, 'k-.', 'LineWidth', 1.5); hold on;
plot(f/1e3, P_clutter, 'g--', 'LineWidth', 1.5); hold on;
plot(f/1e3, P_noise, 'b-.', 'LineWidth', 0.5); 
lgd = legend('目标响应|H(f)|^2', '混响PSD', '噪声PSD');
lgd.FontSize = 14;                          % 图例字号
lgd.Box      = 'off';                       % 去掉边框（可选）
xlabel('频率 kHz', 'FontSize', 18);
ylabel('幅度', 'FontSize', 18);
set(gca, 'FontSize', 17);       % 坐标刻度字号
grid on;


upsample_factor = 10;  % 上采样倍数
N_original = length(f);
N_upsampled = (N_original - 1) * upsample_factor + 1;
% 创建新的频率向量
f_upsampled = linspace(f(1), f(end), N_upsampled);
% 对各个谱进行插值
Pc = interp1(f, P_clutter, f_upsampled, 'pchip');
Pn = interp1(f, P_noise, f_upsampled, 'pchip');
H2 = interp1(f, P_target, f_upsampled, 'pchip');
% 处理可能的NaN值（边界情况）
% Pc(isnan(Pc)) = 0;
% Pn(isnan(Pn)) = 0;
% H2(isnan(H2)) = 0;
f = f_upsampled;


figure;
subplot(2,1,1)
plot(f/1e3, H2, 'k-.', 'LineWidth', 1.5); hold on;
plot(f/1e3, Pc, 'g--', 'LineWidth', 1.5); hold on;
plot(f/1e3, Pn, 'b-.', 'LineWidth', 0.5); 
lgd = legend('目标响应|H(f)|^2', '混响PSD', '噪声PSD');
lgd.FontSize = 14;                          % 图例字号
lgd.Box      = 'off';                       % 去掉边框（可选）
xlabel('频率 kHz', 'FontSize', 18);
ylabel('幅度', 'FontSize', 18);
set(gca, 'FontSize', 17);       % 坐标刻度字号
grid on;
subplot(2,1,2)
plot(f/1e3, 10*log10(H2), 'k-.', 'LineWidth', 1.5); hold on;
plot(f/1e3, 10*log10(Pc), 'g--', 'LineWidth', 1.5); hold on;
plot(f/1e3, 10*log10(Pn), 'b-.', 'LineWidth', 0.5); 
lgd = legend('目标响应|H(f)|^2', '混响PSD', '噪声PSD');
lgd.FontSize = 14;                          % 图例字号
lgd.Box      = 'off';                       % 去掉边框（可选）
xlabel('频率 kHz', 'FontSize', 18);
ylabel('幅度 dB', 'FontSize', 18);
set(gca, 'FontSize', 17);       % 坐标刻度字号
grid on;

% -------------- Test Start ----------------------
% 补0以满足后续波形优化归一化f的条件
% f = linspace(-f(end), f(end), 2*N_upsampled-1);
% H2 = [zeros(1,N_upsampled-1), H2];
% Pc = [zeros(1,N_upsampled-1), Pc];
% Pn = [zeros(1,N_upsampled-1), Pn];

figure;
subplot(2,1,1)
plot(f/1e3, H2, 'k-.', 'LineWidth', 1.5); hold on;
plot(f/1e3, Pc, 'g--', 'LineWidth', 1.5); hold on;
plot(f/1e3, Pn, 'b-.', 'LineWidth', 0.5); 
lgd = legend('目标响应|H(f)|^2', '混响PSD', '噪声PSD');
lgd.FontSize = 14;                          % 图例字号
lgd.Box      = 'off';                       % 去掉边框（可选）
xlabel('频率 kHz', 'FontSize', 18);
ylabel('幅度', 'FontSize', 18);
set(gca, 'FontSize', 17);       % 坐标刻度字号
grid on;
subplot(2,1,2)
plot(f/1e3, 10*log10(H2), 'k-.', 'LineWidth', 1.5); hold on;
plot(f/1e3, 10*log10(Pc), 'g--', 'LineWidth', 1.5); hold on;
plot(f/1e3, 10*log10(Pn), 'b-.', 'LineWidth', 0.5); 
lgd = legend('目标响应|H(f)|^2', '混响PSD', '噪声PSD');
lgd.FontSize = 14;                          % 图例字号
lgd.Box      = 'off';                       % 去掉边框（可选）
xlabel('频率 kHz', 'FontSize', 18);
ylabel('幅度 dB', 'FontSize', 18);
set(gca, 'FontSize', 17);       % 坐标刻度字号
grid on;
% -------------- Test End ----------------------

Be = 4.7;                  % Equivalent Bandwidth Constraint
% Ex = length(f)*2;                  % Signal Energy
Ex = 1e10;
Q_ub = 1;                   % Reverbration Constraint, 实际上限是Q_ub * Ex^2

% target info
v_target = 3;               % m/s
fc = (f_highbound + f_lowbound) / 2;
eta = 1 + 2*v_target / c_v;
f_doppler_target = fc * (eta - 1);
fd_num = round(f_doppler_target / df);

% CM sequence generation params
max_iter = 5000;
tolerance = 1e-5 * Np;
fs_tx_sig = 1e6;
len_gen_sig = 1000;
len_tx_sig = fs_tx_sig * Tp;
sig_Upsample = len_tx_sig / len_gen_sig;
fs_gen_sig = fs_tx_sig / sig_Upsample;

% methods_to_run = {'Flat', 'Waterfill', 'Waterfill_BW', 'Waterfill_Q'};
methods_to_run = {'Waterfill'}; 
results = struct();
for i = 1:length(methods_to_run)
    method = methods_to_run{i};
    fprintf('**** %s 算法 ****\n', method);

    switch method
        case 'Waterfill'
            X_ESD = method_waterfill(H2, Pc, Pn, f, Ex);
        case 'bangbang'
            X_ESD = method_bangbang(H2, Pc, f, Ex);
        case 'Waterfill_BW'
            X_ESD = method_joint(H2, Pc, f, Ex, Be);
        case 'Waterfill_Q'
            X_ESD = method_Qfunc(H2, Pc, Pn, f, Ex, Be, Q_ub, fd_num);
        case 'Waterfill_Comb'
            X_ESD = method_Qfunc_comb(H2, Pc, Pn, f, Ex, fd_nomilized, Be);
        case 'Flat'
            X_ESD = method_constant(f, Ex, 0.95);      % !!!!此处的等效带宽改为1了，SINR更低但是旁瓣表现更好
        otherwise
            warning('未知方法: %s', method);
            continue;
    end
    results.(method).X_ESD = X_ESD;
    results.(method).SINR = SINR_compute(X_ESD, H2, Pc, Pn);
    fprintf('SINR(dB): %.4f\n', 10 * log10(results.(method).SINR));
    fprintf('fd处Q函数值为： %.4f \n',  Qfunc(X_ESD, fd_num, f, Ex));
    fprintf('算法等效带宽： %.4f \n', equivalent_bandwidth(X_ESD, f));
    fprintf('理论能量： %.4f \n', energy(X_ESD, f))
    
    if strcmp(method,'LFM')
        ESD_synthesized = X_ESD;
    else
        [signal, ESD_synthesized] = synthesize_signal_from_ESD(X_ESD.', len_gen_sig, max_iter, tolerance);
    end
    % -------------- Test Start ----------------------
    figure;plot(real(signal));title('Generated Transmit Signal Time Domain');
    figure;plot(f/1e3, ESD_synthesized, 'r-', 'LineWidth', 0.5);title('ESD synthesized');grid on;

    N_gene_sig = length(signal);
    gene_sig_fft = fft(signal);
    gene_sig_fft_shift = fftshift(gene_sig_fft);
    % 计算频率轴
    f_gene_sig = (0:N_gene_sig-1)*(fs_gen_sig/N_gene_sig);
    f_gene_sig_shift = (-N_gene_sig/2:N_gene_sig/2-1)*(fs_gen_sig/N_gene_sig);
    % 计算能量谱密度
    ESD = abs(gene_sig_fft).^2;
    ESD_shift = abs(gene_sig_fft_shift).^2;
    % 绘制能量谱密度
    % figure;
    % plot(f_gene_sig(1:N_gene_sig/2) / 1e3, ESD(1:N_gene_sig/2));
    % xlabel('频率 (kHz)');
    % ylabel('能量谱密度');
    % title('原始生成时域信号的单边能量谱密度');
    figure;
    plot(f_gene_sig_shift / 1e3, ESD_shift);
    xlabel('频率 (kHz)');
    ylabel('能量谱密度');
    title('原始生成时域信号的双边能量谱密度');grid on;
    % -------------- Test End ----------------------

    % 对signal进行上采样
    signal = resample(signal, sig_Upsample, 1).';     % N * 1

    % -------------- Test Start ----------------------
    % figure;plot(real(signal));title('Generated Transmit Signal UPsample Time Domain');
    % figure;plot(f/1e3, ESD_synthesized, 'r-', 'LineWidth', 0.5);title('Generated Transmit Signal UPsample ESD synthesized');
    N_tx_sig = length(signal);
    tx_sig_fft = fft(signal);
    tx_sig_fft_shift = fftshift(tx_sig_fft);
    % 计算频率轴
    f_tx_sig = (0:N_tx_sig-1)*(fs_tx_sig/N_tx_sig);
    f_tx_sig_shift = (-N_tx_sig/2:N_tx_sig/2-1) * (fs_tx_sig/N_tx_sig);
    % 计算能量谱密度
    ESD_shift = abs(tx_sig_fft_shift).^2;
    % 绘制能量谱密度
    figure;
    plot(f_tx_sig_shift / 1e3, ESD_shift);
    xlabel('频率 (kHz)');
    ylabel('能量谱密度');
    title('生成信号上采样的双边能量谱密度');grid on;

    tx_sig_real_fft = fft(real(signal));
    % 计算频率轴
    f_tx_sig = (0:N_tx_sig/2)*(fs_tx_sig/N_tx_sig);
    % 计算能量谱密度
    ESD = abs(tx_sig_real_fft).^2;
    ESD_single_side = ESD(1:N_tx_sig/2+1);
    ESD_single_side(2:end-1) = 2 * ESD_single_side(2:end-1);
    % 绘制能量谱密度
    figure;
    plot(f_tx_sig / 1e3, ESD_single_side);
    xlabel('频率 (kHz)');
    ylabel('能量谱密度');
    title('生成实信号上采样的单边能量谱密度');grid on;
    % -------------- Test End ----------------------

    fprintf('时域波形能量： %.4f \n', sum(abs(signal).^2))
    fprintf('合成ESD能量： %.4f \n', energy(ESD_synthesized, f))
    results.(method).ESD_synthesized = ESD_synthesized;
    results.(method).signal = signal;
    results.(method).autocorr = xcorr(signal, 'normalized'); % 归一化自相关
    results.(method).lags = -length(signal)+1 : length(signal)-1; % 时延序列

    % 可视化未归一化的理论ESD和未归一化的时域波形ESD
    % figure;
    % plot(f, results.(method).X_ESD, 'b--', 'LineWidth', 1.5);
    % hold on;
    % plot(f, results.(method).ESD_synthesized, 'r-', 'LineWidth', 0.5);
    % xlabel('频率');
    % ylabel('幅度');
    % legend('目标ESD', '合成ESD', 'Box', 'off');
    % grid on;

    % 生成信号插值和保存
    signal_real = real(signal);
    % figure;plot(signal_real)
    filename = sprintf('%s_fs1000de3_%dms_fc%dk_B%dk.txt', method, Tp*1e3, fc/1e3, Bw/1e3);
    save(filename, 'signal_real', '-ascii');


end

% ---------------------- 可视化 -------------------------------
colors = {'b','r','m','c','g','y'};  % 自动配色用
% 各自方法波形的ESD
% figure;
% for i = 1:length(methods_to_run)
%     method = methods_to_run{i};
%     plot(f/1e3, results.(method).X_ESD, 'LineWidth', 0.5, 'Color', colors{i}); hold on;
% end
% legend(methods_to_run, 'Location', 'bestoutside');
% title('Waveform Power Spectrum |X(f)|^2');
% xlabel('Frequency kHz');
% ylabel('ESD |X(f)|^2');
% grid on;

%% ------------------------  封装函数    ------------------------------------%


function X_ESD = method_waterfill(H2, Pc, Pn, f, Ex)
% Water filling (max eqution 1)
% Refer to:
% Theory and application of SNR and mutual information matched illumination waveforms
    Bf = sqrt(H2 .* Pn) ./ Pc;
    Df = sqrt(Pn ./ H2);
    % 注水算法：通过二分法确定常数 A = 1 / sqrt(lambda)
    A_low = min(Df); % A 的下界
    A_high = max(Df + Ex ./ Bf); % A 的上界（估计值）
    tolerance = 1e-4; % 容差
    max_iter = 3000; % 最大迭代次数

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

function [ambg, tau, fd] = compute_ambiguity(signal, Fs, PRF)
    % 计算波形的模糊函数
    % 输入:
    %   signal: 时域信号 (列向量)
    %   Fs: 采样率 (Hz)
    %   PRF: 脉冲重复频率 (Hz), 用于多普勒范围归一化
    % 输出:
    %   ambg: 模糊函数矩阵 (多普勒×时延)
    %   tau: 时延向量 (秒)
    %   fd: 多普勒频移向量 (Hz)
    
    N = length(signal);
    tau_max = N/Fs; % 最大时延
    tau = linspace(-tau_max, tau_max, 2*N-1);
    
    % 多普勒范围设置 (±PRF/2)
    fd_max = Fs/25;
    fd_num = 101; % 多普勒点数
    fd = linspace(-fd_max, fd_max, fd_num);
    
    % 计算模糊函数
    ambg = zeros(fd_num, 2*N-1); % 注意这里行列顺序调换
    for k = 1:fd_num
        doppler_signal = signal .* exp(1i*2*pi*fd(k)*(0:N-1)'/Fs);
        ambg(k,:) = abs(xcorr(signal, doppler_signal, 'normalized')).';
    end
end

function [start_idx, end_idx, first_fall_idx] = ...
    findAnalysisInterval(chanvals, threshold, N_analyse, varargin)
%FINDANALYSISINTERVAL  根据阈值寻找第一个“下降-再上升”区间
%
%   [start_idx, end_idx, first_fall_idx] = ...
%       findAnalysisInterval(chanvals, threshold, N_analyse)
%
% 输入
%   chanvals      : 向量或矩阵。若输入矩阵，默认使用第 2 列（与旧脚本一致）。
%   threshold     : 阈值标量
%   N_analyse     : 需要截取的采样点数
%
% 可选名称-值对
%   'Column'      : 指定使用哪一列（避免硬编码第 2 列）
%   'Offset'      : 在找到的 start_idx 上额外平移的量（可正可负）
%
% 输出
%   start_idx     : 第一个上升沿位置
%   end_idx       : start_idx + N_analyse - 1
%   first_fall_idx: 第一个下降沿位置（调试用，可缺省接收）

% ---------------- 解析输入 ----------------
p = inputParser;
addRequired(p,'chanvals');
addRequired(p,'threshold');
addRequired(p,'N_analyse');
addParameter(p,'Column',2,@isnumeric);
addParameter(p,'Offset',0,@isnumeric);
parse(p,chanvals,threshold,N_analyse,varargin{:});
opt = p.Results;

% 若输入向量，强制转为列
if isvector(opt.chanvals)
    opt.chanvals = opt.chanvals(:);
end
col = opt.Column;
v   = opt.chanvals(:,col);

% ---------------- 主逻辑 ------------------
above_threshold = v > opt.threshold;
transitions     = diff(above_threshold);

% 第一个下降沿
first_fall_idx = find(transitions < 0, 1, 'first');
if isempty(first_fall_idx)
    start_idx = [];
    end_idx   = [];
    return
end

% 在该下降沿之后找第一个上升沿
rise_after_fall = find(transitions(first_fall_idx:end) > 0, 1, 'first');
if isempty(rise_after_fall)
    start_idx = [];
    end_idx   = [];
    return
end

start_idx = rise_after_fall + first_fall_idx - 1 + opt.Offset;
start_idx = max(start_idx, 1);                       % 防越界
end_idx   = start_idx + N_analyse - 1;
end
