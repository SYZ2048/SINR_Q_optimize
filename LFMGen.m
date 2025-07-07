clc;
clear all;
fs=1e6; % ²ÉÑùÂÊ
fmin=20e3;
fmax=30e3;
T1=0.5e-3;    % Âö¿í
det=1/fs;
t1=0:det:T1;
W=fmax-fmin;
y=sin(2*pi*fmin*t1 + pi*W/T1*t1.^2);
plot(t1,y);
y=y';
xlswrite('E:\signal\lfm20k_30k_0p5ms.xls',y);