% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��С�������
% ʱ�䣺20180106
% ���������ű�����
% change log�����봦��Ǿ���ȱ���ź�
% (Caution) 
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
% theta = (phi'*phi)*phi'*sig_ls;
% theta = phi\sig_ls;
theta = phi'*sig_ls;

%% �ع��ź�
resample_rete = 100; % [�ɵ�] 
t_full = min(t):1/fs/resample_rete:max(t);
freq_range_re = [fs/N/resample_rete:fs/N/resample_rete:fs/2/resample_rete];
phi = zeros(length(freq_vector)*2+1, length(t_full));
phi(end,:)=1;
for i = 1:length(freq_vector)
    phi(2*i-1,:) = sin(2*pi*freq_vector(i)*t_full);
    phi(2*i, :)  = cos(2*pi*freq_vector(i)*t_full);
end
sig_re = phi'*theta;
figure,plot(t_full, sig_re/max(sig_re))
hold on
plot(t,sig,'g.')
hold off
title('LS��� δ��Լ��')

%%
ft = abs(fftshift(fft(sig)));
L = zeros(N/2+1,1);
L(N/2+1) = theta(end);
for n = 1:N/2
    L(n) = sqrt(theta(2*n-1)^2 + theta(2*n)^2);
end

ft = ft/max(ft);
L = L/max(L);
figure,plot(L),hold on
plot(ft(length(ft)/2+2:length(ft))),hold off


