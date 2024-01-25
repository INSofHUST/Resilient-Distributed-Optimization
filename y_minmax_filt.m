function z_ytemp = y_minmax_filt(i,F,L,ytemp)

% 获得节点i的辅助点
y_i = ytemp(:,i);
z_ytemp = zeros(size(ytemp,1),1);

% 获得节点i邻居的辅助点
n = size(ytemp,2);
neighbor = select_neighbor(L(i,:));
neighbor_index = find(ismember(1:n, neighbor));
neighbor_yi = ytemp(:,neighbor_index);

% 初始化待删集合
[num_dimensions, num_nodes] = size(neighbor_yi);
V_remove = ones(num_dimensions, num_nodes);

% 在每个维度上进行迭代
for dim = 1:num_dimensions
    % 获得该维度上的所有值
    values_dim = neighbor_yi(dim, :);

    % 找到除了节点i本身值的 F-largest and F-smallest 个值的索引
    idx_largest = find(ismember(values_dim,maxk(values_dim(values_dim ~= y_i(dim)), F)));
    idx_smallest = find(ismember(values_dim,mink(values_dim(values_dim ~= y_i(dim)), F)));

%     % 提取可能要移除的值的索引，若值重复，只选一个
%     [~, unique_idx_largest] = unique(values_dim(idx_largest));
%     idx_largest = idx_largest(unique_idx_largest);
% 
%     [~, unique_idx_smallest] = unique(values_dim(idx_smallest));
%     idx_smallest = idx_smallest(unique_idx_smallest);

    % 与xi的值作比较，确定要移除的点
    idx_largest = idx_largest(values_dim(idx_largest) > y_i(dim));
    idx_smallest = idx_smallest(values_dim(idx_smallest) < y_i(dim));
    index_union = union(idx_largest,idx_smallest);
    V_remove(dim,index_union) = 0;

    % 找到保留列向量的序号
    remain_index = find(all(V_remove(dim,:) == 1, 1));

    % 计算剩余状态
    y_mm_i = values_dim(:,remain_index);

    % 自身节点+可信邻居
    y_remain_temp = [y_i(dim) y_mm_i];

    % 动态权重更新
    numneighbor = size(y_mm_i,2);
    numtotal = size(y_remain_temp,2);
    neighbor_weight = ones(1,numneighbor) ./ (numtotal+1);
    total_weight = [1-sum(neighbor_weight) neighbor_weight];

    % 更新辅助点
    z_ytemp(dim,1) = total_weight * y_remain_temp'; 

end
end

