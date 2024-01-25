function [x_dist_i] = dist_filt(i,F,L,ytemp,xtemp)

% 计算接受到的状态到辅助点间的距离
n = size(xtemp,2);
D_ij = sum((xtemp - ytemp(:,i)).^2,1);

% 获得节点i的邻居，及其距离
neighbor = select_neighbor(L(i,:));
neighbor_index = find(ismember(1:n, neighbor));
neighbor_distance = D_ij(neighbor_index);

% 降序排序
[~, sorted_index] = sort(neighbor_distance,'descend');

% 提取最大的F个距离，与Dii比较，移除比Dii大的异常点
D_ii = D_ij(i);
to_remove_index = sorted_index(1:min(F, length(sorted_index)));
to_remove_index = to_remove_index(neighbor_distance(to_remove_index) > D_ii);

% 提取剩余状态
remaining_index = setdiff(1:length(neighbor_index), to_remove_index);
x_dist_i = xtemp(:,neighbor_index(remaining_index));

end

