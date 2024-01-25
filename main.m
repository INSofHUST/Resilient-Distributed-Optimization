% 参考文献：On the Geometric Convergence of Byzantine-Resilient Distributed Optimization Algorithms
% 作者：李云龙
% 时间：2023.11.11
clc;
clear all;

%% 1.初始化节点状态 (i,j,k) 表示第k个节点 j时刻 第i个维度的值
% 分别记录6个常规节点在不同时刻下的演化曲线
x(:,1,1) = [-1 1]';
x(:,1,2) = [-2 2]';
x(:,1,3) = [-3 3]';
x(:,1,4) = [-4 4]';
x(:,1,5) = [-5 5]';
x(:,1,6) = [-6 6]';
x(:,1,7) = [-7 7]';
x(:,1,8) = [-8 8]';
x(:,1,9) = [-9 9]';
x(:,1,10) = [-10 10]';
x(:,1,11) = [-11 11]';
x(:,1,12) = [-12 12]';
x(:,1,13) = [-13 13]';
x(:,1,14) = [-14 14]';
x(:,1,15) = [-15 15]';

% 辅助点
y(:,1,1) = [-1 1]';
y(:,1,2) = [-2 2]';
y(:,1,3) = [-3 3]';
y(:,1,4) = [-4 4]';
y(:,1,5) = [-5 5]';
y(:,1,6) = [-6 6]';
y(:,1,7) = [-7 7]';
y(:,1,8) = [-8 8]';
y(:,1,9) = [-9 9]';
y(:,1,10) = [-10 10]';
y(:,1,11) = [-11 11]';
y(:,1,12) = [-12 12]';
y(:,1,13) = [-13 13]';
y(:,1,14) = [-14 14]';
y(:,1,15) = [-15 15]';

%% 2.描述图的拓扑结构
% e.g 全连接矩阵
[d,~,n] = size(x);
% 邻接矩阵，只需要根据图拓扑修改A即可
A = AdjacencyMatrix(n);
% 度矩阵D 
D = diag(sum(A,2)); 
% 拉普拉斯矩阵 
L = D - A; 
% 行随机权重矩阵，满足行和为1即可，参照Nedic 2009 TAC 
W = Weight(L); 

% 对抗节点集合 
Adversary = [14 15]; 
% 常规节点集合 
Regularnodes = setdiff(1:n, Adversary);
numReg = length(Regularnodes);

%% 3.变量储存
for i = 1:n
    xtemp(:,i) = x(:,1,i);      % 中间变量，记录第k时刻，所有结点的状态
    xsyn(:,i) = x(:,1,i);       % 中间变量，记录第k+1时刻，所有节点同步更新后的状态

    ytemp(:,i) = y(:,1,i);      % 辅助点变量，记录第k时刻，所有结点的状态
    ysyn(:,i) = y(:,1,i);       % 辅助点变量，记录第k+1时刻，所有节点同步更新后的状态
end

%% 4.仿真迭代
T = 5000;                % 设置迭代次数
stepC = 1;                % 设置迭代倍数
F = length(Adversary);    % 对抗节点数目上界

for k = 1:T
    
    % 迭代步长 衰减的
    stepK = 1/(k+1);         

    for i = 1:n                          % 并行计算，每个结点分别同时计算
        
        % 节点i的状态值
        x_i = xtemp(:,i);

%         % 算法1：不加工
%         x_mm_i = x_minmax_filt(F,xtemp,x_i);
        
        % 算法2：K.K
          x_mm_i = SDMMFD(i,F,L,ytemp,xtemp);
        
        
        % 动态权重更新
        x_remain_temp = [x_i x_mm_i];
        numneighbor = size(x_mm_i,2);
        numtotal = size(x_remain_temp,2);
        neighbor_weight = ones(1,numneighbor) ./ (numtotal);
        total_weight = [1-sum(neighbor_weight) neighbor_weight];

        % 加权
        z_xtemp = total_weight * x_remain_temp';          
        z_xtemp = z_xtemp';

        % 节点i 下一步到达的状态
        xsyn(:,i) = z_xtemp - stepK * grad(i,z_xtemp); 
        
        % 节点i 辅助点下一时刻的状态
        z_ytemp = y_minmax_filt(i,F,L,ytemp);
        ysyn(:,i) = z_ytemp;

    end

    % 攻击者发起攻击
    xsyn(:,Adversary(1)) = 2*k+10*rand(1,2)';
    ysyn(:,Adversary(1)) = 5*rand(1,2)';

    xsyn(:,Adversary(2)) = 10*rand(1,2)';
    ysyn(:,Adversary(2)) = 10*rand(1,2)';

    % 读取状态
    for i = 1:n
        xtemp(:,i) = xsyn(:,i);
        x(:,k+1,i) = xsyn(:,i);
        ytemp(:,i) = ysyn(:,i);
        y(:,k+1,i) = ysyn(:,i);
    end

end

%% 5.计算真实值
% 定义目标函数
fun_i = @(i, x) (x(1) + i)^2 + (x(2) - i)^2;
% 定义目标函数，即所有节点函数的和
fun_sum = @(x) sum(arrayfun(@(i) fun_i(i, x), [Regularnodes]));
% 初始猜测值
x0 = zeros(2,1);
% 求解无约束最小化问题，即求解和函数真实最优解
options = optimset('Display', 'iter'); % 设置显示迭代过程
[x_min, fval] = fminunc(fun_sum, x0, options);
% 输出
disp(x_min)
disp(['极小值 f(x): ', num2str(fval)])

%% 6.一致性分析
% 第m维 一致性收敛曲线
figure
% 设置每行显示的图例个数
r = 4;

for k = 1:d
    subplot(d, 1, k)
    % 使用 for 循环遍历向量
    for i = 1:length(Regularnodes)
        plot(1:T, x(k, 1:T, Regularnodes(i)), '-.', 'linewidth', 1.5, 'DisplayName', ['Agent ', num2str(Regularnodes(i))]);
        hold on;
    end
    plot(1:T, ones(size(1:T)) * x_min(k), 'k.-', 'linewidth', 2, 'DisplayName', ['min value']);
    hold on;
    title(['Variable ', num2str(k), ' Trajectories']);
    xlabel('Time');
    ylabel(['Variable ', num2str(k)]);
    % 显示图例
    hLegend = legend('show');
    % 去除图例的边框
    set(hLegend, 'Box', 'off');
    % 调整图例的位置和布局
    legend('NumColumns', r, 'Location', 'best');
end

%% 7.最优值误差
% 计算不同时刻 fun_sum 的值
fun_sum_values = zeros(1, T);
for t = 1:T
    % 提取当前时刻的节点状态
    x_t = squeeze(x(:, t, Regularnodes));
    % 计算当前时刻的 fun_sum
    fun_sum_t = sum(arrayfun(@(i) fun_i(i, x_t), Regularnodes));
    % 存储结果
    fun_sum_values(t) = fun_sum_t;
end

%% 8.准确性
constant_value = fun_sum(x_min);
accuracy_value = (1 - (fun_sum_values - constant_value)/constant_value)*100;

%% 9.最优误差
error = fun_sum_values - constant_value;
figure;
% 选择以5为倍数的时刻
selected_times = 1:T;
plot(selected_times, error(selected_times), '-', 'LineWidth', 1.5, 'DisplayName', 'Ours');
hold on;

% 修改图表属性
title('The error of optimal solution', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Number of iterations $K$', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

% 设置 y 轴刻度
set(gca, 'YScale', 'log');
% 添加图例，并设置位置和字体大小
legend('Orientation', 'vertical', 'Location', 'northeast', 'FontSize', 12, 'Interpreter', 'latex', 'Box', 'off', 'Position', [0.7, 0.8, 0.1, 0.1]);