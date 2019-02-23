% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最小二乘替代
% 时间：20180106
% 附属函数脚本：无
% change log：对三种求解方法的分析对比
% (Caution) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
fs = 512;
% 均匀时间
t = 1/fs:1/fs:2;
N = length(t);
sig = sin(2*pi*(20*t'.^2 + 10*t'));
%% 制造缺失区间
jump_op = 319; % [可调] 
jimp_ed = 320; % [可调] 
t_ls = t([1:jump_op,jimp_ed:end]);
sig_ls = sig([1:jump_op,jimp_ed:end]);
N_ls = length(t_ls);
continue_rate = 4/4; % [可调] 
%% 无约束最小二乘数值解
% 问题：直流分量出现在第一行会导致结果出错
% 求矩阵逆时报错并产生不理想的结果，论文中基函数顺序建议调整
freq_range = [fs/N:fs/N:fs/2];
freq_vector = freq_range(1:floor(length(freq_range)*continue_rate));
phi = ones(N_ls, length(freq_vector)*2+1);
for n = 1:length(freq_vector)
    phi(:,2*n-1) = sin(2*pi*freq_vector(n)*t_ls);
    phi(:,2*n) = cos(2*pi*freq_vector(n)*t_ls);
end
% 求解LS方法
% 问题：使用求解方法1，方法3效果和FFT缺失信号置0效果相同，具体表现为恢复信号缺失部分几乎为0
% 方法2是左除，结果与张心亮学长程序相同，出现缺失部分过拟合问题
% inv(phi'*phi)的结果是仅有对角线有数值的单位矩阵
base1 = (phi'*phi)*phi';
base2 = inv(phi'*phi)*phi';
base2 = base2./max(base2,[],2);
theta1 = base1*sig_ls; 
theta2 = phi\sig_ls;
theta3 = phi'*sig_ls;
%% 重构信号
resample_rete = 1; % [可调] 
t_full = min(t):1/fs/resample_rete:max(t);
freq_range_re = [fs/N/resample_rete:fs/N/resample_rete:fs/2/resample_rete];
phi = zeros(length(freq_vector)*2+1, length(t_full));
for i = 1:length(freq_vector)
    phi(2*i-1,:) = sin(2*pi*freq_vector(i)*t_full);
    phi(2*i, :)  = cos(2*pi*freq_vector(i)*t_full);
end
sig_re1 = phi'*theta1;
sig_re2 = phi'*theta2;
sig_re3 = phi'*theta3;
%% 重构信号与原信号时域对比
figure,plot(t_full, sig_re1/max(sig_re1))
hold on
plot(t,sig,'g.')
hold off
title('LS未加约束')
%% 重构信号与原信号频域对比
ft = abs(fftshift(fft(sig)));
L1 = zeros(N/2+1,1);
L2 = zeros(N/2+1,1);
L3 = zeros(N/2+1,1);
L1(N/2+1) = theta1(end);
L2(N/2+1) = theta2(end);
L3(N/2+1) = theta3(end);
for n = 1:N/2
    L1(n) = sqrt(theta1(2*n-1)^2 + theta1(2*n)^2);
    L2(n) = sqrt(theta2(2*n-1)^2 + theta2(2*n)^2);
    L3(n) = sqrt(theta3(2*n-1)^2 + theta3(2*n)^2);
end
% 归一化
ft = ft/max(ft);
L1 = L1/max(L1);
L2 = L2/max(L2);
L3 = L3/max(L3);
figure,plot(L1),hold on
plot(L2),hold on
plot(L3),hold on
plot(ft(length(ft)/2+2:length(ft))),hold off
legend('base','左除','转置', 'FFT')

