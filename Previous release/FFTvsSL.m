% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��С�������FFT
% �о�ȱʧ���������LS�����쳣�Ĵ���취
% ���������ű���
%           FFTvsSL_ParsevalConstraint.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;%close all;

e = 2.7183;
fs = 512;
t = [1/fs:1/fs:2.0];

%ʵ���ź�
h = sin(2*pi*(20*t.^2 + 10*t))'; % ���Ե�Ƶ nN=1.3 mM=2.3 sh=560
N = length(h);


freq_range = [fs/N:fs/N:fs/2];
freq_range_FFT = freq_range;
H = abs((fft(h)));
H = H/length(h);



%% ����ȱʧ����
jump_op = 280; % [�ɵ�] 
jimp_ed = 320; % [�ɵ�] 
t_ls = t([1:jump_op,jimp_ed:end]);
h_ls = h([1:jump_op,jimp_ed:end]);
N_ls = length(t_ls);

freq_range_LS = [fs/N:fs/N:fs/2];
continue_rate = 4/4; % [�ɵ�] 

%% ���Լ��
do_ParsevalConstraint = 1;
do_blast = 0;
do_NewtonLagrange = 1;
do_SQP = 0;
param_y_div_x = 2;
%% ��Լ����С������ֵ��
freq_vector = freq_range(1:floor(length(freq_range)*continue_rate));
phi_mat = zeros(N_ls, length(freq_vector)*2+1);
phi_mat(:,end)=1;
for iii = 1:length(freq_vector)
    phi_mat(:,2*iii-1) = sin(2*pi*freq_vector(iii)*t_ls);
    phi_mat(:,2*iii)  = cos(2*pi*freq_vector(iii)*t_ls);
end
LS_temp = pinv(phi_mat)*h_ls;


%% ��Լ����������
if sum(LS_temp.^2)/2 > 1.1*sum(h_ls.^2)/length(h_ls)
    disp('LS��������')
end


if do_ParsevalConstraint == 1
    disp('���������')
%     load LS_param_last
    if do_blast == 1
        FFTvsSL_PC_blast;
    elseif do_NewtonLagrange == 1
        FFTvsSL_PC_NewtonLagrange;
    elseif do_SQP == 1
        Copy_of_FFTvsSL_PC_SQP;
    end
end


%%

% load LS_param_last
% LS_temp = LS_param_last

LS = zeros(length(freq_vector)+1,1);
if do_ParsevalConstraint == 1
    LS(1) = LS_temp_new(end);
    for ii = 1:length(freq_vector)
        LS(ii+1) = sqrt(LS_temp_new(2*ii-1).^2 + LS_temp_new(2*ii).^2);
    end
else
    LS(1) = LS_temp(end);
    for ii = 1:length(freq_vector)
        LS(ii+1) = sqrt(LS_temp(2*ii-1).^2 + LS_temp(2*ii).^2);
    end
end

% �ع�
resample_rete = 1; % [�ɵ�] 
t_full = min(t):1/fs/resample_rete:max(t);
freq_range_re = [fs/N/resample_rete:fs/N/resample_rete:fs/2/resample_rete];
Phi_mat = zeros(length(freq_vector)*2+1, length(t_full));
Phi_mat(end,:)=1;
for iii = 1:length(freq_vector)
    Phi_mat(2*iii-1,:) = sin(2*pi*freq_vector(iii)*t_full);
    Phi_mat(2*iii, :)  = cos(2*pi*freq_vector(iii)*t_full);
end
if do_ParsevalConstraint == 1
    h_re = Phi_mat' * LS_temp_new;
else
    h_re = Phi_mat' * LS_temp;
end

% h_re(find((abs(h_re)-1.2)>0)) = 0;
H_re = abs(fft(h_re));
H_re = H_re/length(h_ls);
% subplot(1,2,1)

if do_ParsevalConstraint == 1
    figure(3)
    plot([0, freq_range], H(1:length(freq_range)+1)*2)
    xlabel('frequency/Hz')
    hold on
    plot([0,freq_vector], LS(1:end), 'c.')
    xlabel('frequency/Hz')
    % plot([0, freq_range], H_re(1:length(freq_range)+1),'m.')
    legend('FFT', 'LS','FFT rebuild')
    plot([0,freq_vector], LS(1:end), 'c')
    % plot([0, freq_range], H_re(1:length(freq_range)+1),'m--')
    title('FFT vs LS ���Լ��')
    hold off
    
    
    % subplot(1,2,2)
    figure(4)
    plot(t_full, h_re)
    hold on
    plot(t,h,'g.')
    plot(t_ls,h_ls,'r.')
    hold off
    title('LS��� ���Լ��')
    
else
    figure(1)
    plot([0, freq_range], H(1:length(freq_range)+1)*2)
    xlabel('frequency/Hz')
    hold on
    plot([0,freq_vector], LS(1:end), 'c.')
    xlabel('frequency/Hz')
    % plot([0, freq_range], H_re(1:length(freq_range)+1),'m.')
    legend('FFT', 'LS','FFT rebuild')
    plot([0,freq_vector], LS(1:end), 'c')
    % plot([0, freq_range], H_re(1:length(freq_range)+1),'m--')
    title('FFT vs LS δ��Լ��')
    hold off
    
    
    % subplot(1,2,2)
    figure(2)
    plot(t_full, h_re)
    hold on
    plot(t,h,'g.')
    plot(t_ls,h_ls,'r.')
    hold off
    title('LS��� δ��Լ��')
end



