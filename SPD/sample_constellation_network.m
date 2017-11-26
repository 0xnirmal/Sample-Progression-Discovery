function [adj,adj_mst,avg_edge_len] = sample_constellation_network(data,average_degree)
% each column of data is one sample

dist_matrix = squareform(pdist(data','cityblock'))/size(data,1); dist_matrix = dist_matrix + eye(size(dist_matrix))*max(max(dist_matrix))*2;
[adj,adj2, cost_value] = mst_from_dist_matrix(dist_matrix); % get mst
adj_mst = adj;
num_edges_needed = average_degree*size(adj,1)/2-sum(sum(adj))/2; % to achieve average_degree, we need in total average_degree*size(adj,1)/2 edges, we have sum(sum(adj))/2 edges, so the difference is the number we need
if num_edges_needed>=0
    tmp = dist_matrix; 
    tmp(adj==1) = max(max(tmp)); 
    [Y,I] = sort(tmp(:),'ascend');
%     adj(tmp<=Y(num_edges_needed*2))=1;
    adj = adj + (tmp<=Y(num_edges_needed*2));
end
avg_edge_len = full(sum(sum((dist_matrix.*adj)))/(sum(sum(adj))));

return
