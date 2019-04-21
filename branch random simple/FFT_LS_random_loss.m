% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最小二乘与FFT对比试验   非均匀采样！
% 时间：201801019
% 附属函数脚本：无
% change log：
% (Caution) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
fs = 512;
% 均匀时间
t =  [1/fs:1/fs:1.0];
N = length(t);
sig = sin(2*pi*(20*t'.^2 + 10*t'));
%%
%%制造非均匀采样
k = 5;
sig_new = sig;
sig_new1 = sig;
t_ls = t;
missing_point = round(1+512.*rand(k,1));
sig_new(missing_point) = []; 
sig_new1(missing_point) = 0;
t_ls(missing_point)= [];
fft_sig_new1 = abs(fftshift(fft(sig_new1)));
figure,plot(fft_sig_new1);

sig_ls = sig_new;
N_ls = length(t_ls);
continue_rate = 4/4;
%% 无约束最小二乘数值解
freq_range = [fs/N_ls:fs/N_ls:fs/2];
freq_vector = freq_range(1:floor(length(freq_range)*continue_rate));
phi = ones(N_ls, length(freq_vector)+1);
for n = 1:length(freq_vector)
    phi(:,n) = exp(1i*2*pi*freq_vector(n)*t_ls);
end
% theta = (phi'*phi)*phi'*sig_ls;
% theta = phi\sig_ls;
theta = phi'*sig_new;

%% 重构信号
resample_rete = 1; % [可调] 
t_full = min(t):1/fs/resample_rete:max(t);
freq_range_re = [fs/N/resample_rete:fs/N/resample_rete:fs/2/resample_rete];
phi = ones(length(t_full),length(freq_vector)+1);
for i = 1:length(freq_vector)
    phi(:,i) = exp(1i*2*pi*freq_vector(i)*t_full);
end
sig_re = phi*theta;
real_sig_re = real(sig_re)/max(real(sig_re));
figure,plot(t_full, real_sig_re)
hold on
plot(t,sig,'g.')
hold off
title('LS拟合 未加约束')
%%
ft = abs(fftshift(fft(sig_new)));
ft = ft/max(ft);
theta = theta/max(theta);
figure,plot(1:length(theta), abs(theta)),hold on
plot(0:1:(length(ft)/2),ft(length(ft)/2:length(ft))),hold off
title('非均匀采样信号的LS 对比 FFT')
%%
ftt = abs(fftshift(fft(sig)));
ftt = ftt/max(ftt);
for i = 1:511
    ftt1(i) = ftt(i+1);
end
figure,plot(1:length(theta), abs(theta)),hold on
plot(0:1:(length(ftt1)/2),ftt1(length(ftt1)/2:length(ftt1))),hold off
title('非均匀采样信号的LS')
%%
ftt = abs(fftshift(fft(sig)));%曲线平滑！！！
ftt = ftt/max(ftt);
ffft = abs(fftshift(fft(sig_new1,512)));
ffft = ffft/max(ffft);
% figure,plot(ffft);hold on 
% plot(ftt);hold off
% title('非均匀采样信号的FFT11')
figure,plot(0:1:(length(ffft)/2), ffft(length(ffft)/2:length(ffft))),hold on
plot(0:1:(length(ftt)/2),ftt(length(ftt)/2:length(ftt))),hold off
title('非均匀采样信号的FFT')