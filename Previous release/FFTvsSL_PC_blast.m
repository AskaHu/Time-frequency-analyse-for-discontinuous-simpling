% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFTvsSL�Ӻ���
% ���������߶�����Լ����С�����Ż����
% �漴���Ʒ�
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ������� OP
        range_parseval = sqrt(sum(h_ls.^2)/length(h_ls)*2);
%         err2_last = 1e10;
%         LS_param_last = LS_temp;
        load LS_param_last
        while 1
            LS_param = rand(size(LS_temp)) + 100*LS_param_last;
            LS_param = LS_param/sqrt(sum(LS_param.^2))*range_parseval;
            
            
            err2 = sum((phi_mat'*LS_param - h_ls).^2);
            
            
            
            
            
            if err2 >= err2_last
                continue;
            else
                err2_last
                err2
                
                
                if err2_last/err2<1.000000001
                    save LS_param_last LS_param_last err2_last
                    break;
                end
                err2_last = err2;
                LS_param_last = LS_param;
            end

            
        end
        
        % �漴���� ED