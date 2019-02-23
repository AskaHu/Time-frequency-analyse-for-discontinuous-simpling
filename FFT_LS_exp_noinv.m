% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��С�������
% ʱ�䣺20180106
% ���������ű�����
% change log����������ⷽ���ķ����Ա�
% (Caution) �������㲻����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
fs = 512;
% ����ʱ��
t = 1/fs:1/fs:2;
N = length(t);
sig = sin(2*pi*(20*t'.^2 + 10*t'));
%% ����ȱʧ����
jump_op = 300; % [�ɵ�] 
jimp_ed = 320; % [�ɵ�] 
t_ls = t([1:jump_op,jimp_ed:end]);
sig_ls = sig([1:jump_op,jimp_ed:end]);
sig_zero = sig;
sig_zero(jump_op+1:jimp_ed-1) = 0;
N_ls = length(t_ls);
continue_rate = 4/4; % [�ɵ�] 
%% ��Լ����С������ֵ��
freq_range = [fs/N:fs/N:fs/2];
freq_vector = freq_range(1:floor(length(freq_range)*continue_rate));
phi = ones(N_ls, length(freq_vector)+1);
for n = 1:length(freq_vector)
    phi(:,n) = exp(1i*2*pi*freq_vector(n)*t_ls);
end
% ���LS����
% Ŀǰֻ�з���3���ã��������������̫�ֲ��������ж�Ϊ����������������
% theta = inv(phi'*phi)*phi'*sig_ls;
% theta = phi\sig_ls;
theta = phi'*sig_ls;
%% �ع��ź�
resample_rete = 100; % [�ɵ�] 
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
plot(t_ls,sig_ls,'r.')
legend('�ع��ź�', 'ԭ�ź�')
hold off
title('EXP LSδ��Լ��')
%% �ع��ź���ԭ�ź�Ƶ��Ա�
ft = abs(fftshift(fft(sig)));
ft_zero = abs(fftshift(fft(sig_zero)));
% ��һ��
ft = ft/max(ft);
theta = theta/max(theta);
ft_zero = ft_zero/max(ft_zero);
figure,plot(1:length(theta), abs(theta)),hold on
plot(0:1:(length(ft)/2),ft(length(ft)/2:length(ft))),hold on
plot(0:1:(length(ft_zero)/2),ft_zero(length(ft_zero)/2:length(ft_zero))),hold off
legend('Our Method', 'FFT δȱʧ','FFT Zero')


