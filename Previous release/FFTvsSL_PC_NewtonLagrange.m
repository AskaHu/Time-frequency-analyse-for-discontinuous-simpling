% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFTvsSL子函数
% 利用帕塞瓦尔定理约束最小二乘优化结果
% 牛顿-拉格朗日法
% 局部二次收敛
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = phi_mat;
y = h_ls;
AtA = A'*A;
Aty = A'*y;
norm_y = norm(y)*sqrt(2/length(y)) * param_y_div_x;

k = 1;
x = LS_temp / norm(LS_temp)*norm_y;
x_len = length(x);
miu = 0.5;


round_max = 400;
z = zeros(x_len+length(miu), round_max);
p = z;
rho = 0.5;
gamma = 0.4;
mk = 0;
tic
while k<round_max
    z(:,k) = [x; miu];
    if k>2
        a(:,k-1) = z(:,k) - z(:,k-1);
    end
    Hesse = [[2*AtA*miu-2*eye(x_len), 2*AtA*x-2*Aty];[(2*AtA*x-2*Aty)', 0]];
    dL = [2*AtA*x*miu-2*Aty*miu-2*x;norm(A*x-y)];
    p(:,k) = Hesse\(-dL);
    
    % Amijio check
    m = floor(mk/100);
    mk = m;
    x_add = p(1:x_len,k);
    miu_add = p(x_len+1:end,k);
    
    while m<=100000
        if norm([2*AtA*(x+x_add*rho.^m)*(miu+miu_add*rho.^m)-2*Aty*(miu+miu_add*rho.^m)-2*(x+x_add*rho.^m);...
                -norm(x+x_add*rho.^m)+norm_y]) < sqrt(1-gamma*rho^m)*norm(dL)
            mk = m;
            break
        end
        if m<50
            m = m+1;
        elseif m<400
            m = m+10;
        else
            m = m+100;
        end
        mk = m;
    end
    if mk >= 100000
        mk
        break
    end
    x = x+rho^mk * p(1:x_len,k);
    miu = miu+rho^mk * p(x_len+1:end,k);
    k = k+1;
end
toc


err2_min_SQP = norm(A*x-y)
norm_LS = norm(x)
LS_temp_new = x;





