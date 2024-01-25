function [x_mm_i] = x_minmax_filt(F,x_dist_i,x_i)

[num_dimensions, num_nodes] = size(x_dist_i);

% 初始化待删集合
V_remove = ones(num_dimensions, num_nodes);

% 在每个维度上进行迭代
for dim = 1:num_dimensions
    % 获得该维度上的所有值
    values_dim = x_dist_i(dim, :);

    % 找到除了节点i本身值的 F-largest and F-smallest 个值
    idx_largest = find(ismember(values_dim,maxk(values_dim(values_dim ~= x_i(dim)), F)));
    idx_smallest = find(ismember(values_dim,mink(values_dim(values_dim ~= x_i(dim)), F)));

    % 与xi的值作比较，确定要移除的点
    idx_largest = idx_largest(values_dim(idx_largest) > x_i(dim));
    idx_smallest = idx_smallest(values_dim(idx_smallest) < x_i(dim));
    index_union = union(idx_largest,idx_smallest);

    V_remove(:,index_union) = 0;
end

% 找到保留列向量的序号
remain_index = find(all(V_remove == 1, 1));

% 计算剩余状态
x_mm_i = x_dist_i(:,remain_index);

end

