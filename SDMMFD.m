function [x_mm_i,z_ytemp] = SDMMFD(i,F,L,ytemp,xtemp)
% K.K 基于辅助点过滤

% 节点i的状态值
x_i = xtemp(:,i);

% x_dist_filter
x_dist_i = dist_filt(i,F,L,ytemp,xtemp);

% x_minmax_filter
x_mm_i = x_minmax_filt(F,x_dist_i,x_i);

end

