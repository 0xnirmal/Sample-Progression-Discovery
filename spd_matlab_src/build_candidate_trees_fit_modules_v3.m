function [candidate_trees, p_value] = build_candidate_trees_fit_modules_v3(data,idx, permutation_iter)
% this is\was part of the file analyze_progression
% when making the GUI, it seems good to break that file into parts
% later, I may use this function to replace the corresponding part in 
% analyze_progression, to make it look cleaner


iter = permutation_iter;
% candidate_trees = zeros(size(data,2),size(data,2),max(idx));
% candidate_graphs = zeros(size(data,2),size(data,2),max(idx));
candidate_trees=[];
for i=1:max(idx)
    candidate_trees{i} = sparse(zeros(size(data,2),size(data,2)));
    candidate_graphs{i} = sparse(zeros(size(data,2),size(data,2)));
end
for i=1:max(idx), 
    display(['constructing tree ', num2str(i)]);
    module_size(i) = sum(idx==i); 
    dist_matrix = squareform(pdist(data(idx==i,:)','cityblock'));
    [adj_graph,adj_mst,avg_edge_len] = sample_constellation_network(data(idx==i,:),3);
    candidate_trees{i} = adj_mst;
    candidate_graphs{i} = adj_graph;
    drawnow;
end

fprintf('\n\n finished building candidate trees for %d modules \n\n', max(idx));

p_value = [];
earth_mover_dist=[];
for j=1:max(idx) % loop over the modules
    dist_matrix = squareform(pdist(data(idx==j,:)','cityblock'));
    num_bins_ = 20;
    [N,X] = hist(squareform(dist_matrix),num_bins_);
    P = N/sum(N);
    c = squareform(pdist((1:num_bins_)','cityblock')); c = c - triu(c);
    for i=1:max(idx) % loop over the trees
        fprintf('fitting candidate trees and modules, pairs %4d, %4d\n',i,j); drawnow;
        adj_mst = full(candidate_trees{i});
        adj_graph = full(candidate_graphs{i});
%         [Q] = hist(dist_matrix(triu(adj_mst,1)==1),X); Q = Q/sum(Q);
%         [dist,F] = EarthMoverDist(P,Q,c);
        [Q] = hist(dist_matrix(triu(adj_graph,1)==1),X); Q = Q/sum(Q);
        [dist,F] = EarthMoverDist(P,Q,c);
        earth_mover_dist(i,j) = dist;
    end
end
earth_mover_dist(earth_mover_dist<=0) = min(earth_mover_dist(earth_mover_dist>0)); % the larger, the better fit
[Y,I] = sort(earth_mover_dist(:),'descend');
dist_rank(I) = 1:length(I);
p_value = reshape(dist_rank/length(dist_rank(I)),size(earth_mover_dist,1),size(earth_mover_dist,2)); 


return

