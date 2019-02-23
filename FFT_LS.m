% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��С�������
% ʱ�䣺20180106
% ���������ű�����
% change log����������ⷽ���ķ����Ա�
% (Caution) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
fs = 512;
% ����ʱ��
t = 1/fs:1/fs:2;
N = length(t);
sig = sin(2*pi*(20*t'.^2 + 10*t'));
%% ����ȱʧ����
jump_op = 319; % [�ɵ�] 
jimp_ed = 320; % [�ɵ�] 
t_ls = t([1:jump_op,jimp_ed:end]);
sig_ls = sig([1:jump_op,jimp_ed:end]);
N_ls = length(t_ls);
continue_rate = 4/4; % [�ɵ�] 
%% ��Լ����С������ֵ��
% ���⣺ֱ�����������ڵ�һ�лᵼ�½������
% �������ʱ��������������Ľ���������л�����˳�������
freq_range = [fs/N:fs/N:fs/2];
freq_vector = freq_range(1:floor(length(freq_range)*continue_rate));
phi = ones(N_ls, length(freq_vector)*2+1);
for n = 1:length(freq_vector)
    phi(:,2*n-1) = sin(2*pi*freq_vector(n)*t_ls);
    phi(:,2*n) = cos(2*pi*freq_vector(n)*t_ls);
end
% ���LS����
% ���⣺ʹ����ⷽ��1������3Ч����FFTȱʧ�ź���0Ч����ͬ���������Ϊ�ָ��ź�ȱʧ���ּ���Ϊ0
% ����2������������������ѧ��������ͬ������ȱʧ���ֹ��������
% inv(phi'*phi)�Ľ���ǽ��жԽ�������ֵ�ĵ�λ����
base1 = (phi'*phi)*phi';
base2 = inv(phi'*phi)*phi';
base2 = base2./max(base2,[],2);
theta1 = base1*sig_ls; 
theta2 = phi\sig_ls;
theta3 = phi'*sig_ls;
%% �ع��ź�
resample_rete = 1; % [�ɵ�] 
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
%% �ع��ź���ԭ�ź�ʱ��Ա�
figure,plot(t_full, sig_re1/max(sig_re1))
hold on
plot(t,sig,'g.')
hold off
title('LSδ��Լ��')
%% �ع��ź���ԭ�ź�Ƶ��Ա�
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
% ��һ��
ft = ft/max(ft);
L1 = L1/max(L1);
L2 = L2/max(L2);
L3 = L3/max(L3);
figure,plot(L1),hold on
plot(L2),hold on
plot(L3),hold on
plot(ft(length(ft)/2+2:length(ft))),hold off
legend('base','���','ת��', 'FFT')

