clc;
clear all;
fc = 100e3;
fs = 3*fc; % ²ÉÑùÂÊ
Bw = 40e3;
% fmin=20e3;
% fmax=30e3;
fmin = fc - Bw/2;
fmax = fc + Bw/2;
T1 = 20e-3;    % Âö¿í
det=1/fs;
t1=0:det:T1;
y=sin(2*pi*fmin*t1 + pi*Bw/T1*t1.^2);
plot(t1,y);
y=y';       % 6001 * 1
xlswrite('C:\DATA\code_MATLAB\SINR_Q_optimize\lfm80k_120k_20ms_fs300e3.xls',y);

%% 
clc
clear ;
close all
%85-115   90-120  90-110
fs = 200e3;     %cai yang lv
T = 200e-3;        %mai kuan
Nw = round(fs*T);
t = (0 : Nw-1) / fs;
B = 10e3;       %dai kuan
k = B / T;

f0_1 = 15e3;       %zhong xin pin lv
fl_1 = f0_1 - 0.5*B;
fh_1 = f0_1 + 0.5*B;

lfm_1 = sin(2*pi*fl_1*t +pi*k*t.^2);
lfm_1(end) = 0;
lfm_1 = lfm_1.';
figure;plot(lfm_1);

lfm_1_fft = fft(lfm_1,fs);
figure;plot(abs(lfm_1_fft));
xlim([f0_1-B,f0_1+B]);
% spectrogram(lfm_1,512,256,fs,fs,'yaxis');
% ylim([60,120])

save('lfm_fs200e3_200ms_fc15k_B10k.txt','lfm_1','-ascii');
figure;
spectrogram(lfm_1,256,128,100,fs,'yaxis');