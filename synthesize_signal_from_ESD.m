function [a_final, ESD_synthesized] = synthesize_signal_from_ESD(ESD_desired, M, max_iter, tolerance)
    % 必须满足：M < N
    % 输入:
    % ESD_desired: 目标ESD (N x 1)
    % M: 时域信号的长度
    % max_iter: 最大迭代次数
    % tolerance: 收敛容忍度 (如0.01*M)
    % 输出:
    % a_final: 合成的时域信号 (M x 1) 
    % ESD_synthesized: 合成信号的ESD (N x 1)

    N = length(ESD_desired);
    S = sqrt(ESD_desired); % 目标幅度谱

    % 初始化相位为0
    phi = zeros(N, 1);
    x = S .* exp(1i * phi); % 初始x向量

    % 创建DFT矩阵
    W = dftmtx(N);
    W = W(:, 1:M); % 取前M列

    a_prev = zeros(M, 1); % 用于存储上一次迭代的结果
    for iter = 1:max_iter
        % 计算最小二乘估计
        a_hat = (1/N) * W' * x;

        % 对a_hat进行归一化，强制幅度为1
        a = exp(1i * angle(a_hat));

        % 计算DFT
        A = W * a;

        % 更新相位
        phi = angle(A);

        % 更新x向量
        x = S .* exp(1i * phi);

        % 检查收敛性
        delta = sum(abs(a - a_prev));
        if delta < tolerance
            fprintf('算法在 %d 次迭代后收敛。\n', iter);
            break;
        end
        a_prev = a;
    end

    a_final = a;
    ESD_synthesized = abs(A).^2;
    ESD_synthesized = ESD_synthesized / max(ESD_synthesized);
end