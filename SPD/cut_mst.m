function [idx] = cut_mst(adj,data,threshold)

% threshold = 0.3;
components = sparse(ones(size(adj,1),1));
if avg_center_gene_corr(data(components(:,1)==1,:))>threshold
    isleaf = 1; 
    isactive=0;
else
    isleaf = 0;
    isactive=1;
end

idx_tmp = zeros(size(adj,1),1);
while sum(isactive==1)~=0
    component_ind = find(isactive==1,1);
    [component_ind, length(isactive)]
    [full(sum(components(:,isleaf==1)))]
%     if  component_ind==9
%         1
%     end
    component_members = find(components(:,component_ind)==1);
    if length(component_members)==1
        isleaf(component_ind) = 1; 
        isactive(component_ind)=0;
        continue;
    end
    if avg_center_gene_corr(data(component_members,:))>threshold
        isleaf(component_ind) = 1; 
        isactive(component_ind)=0;
        continue;
    end
    [idx_tmp] = partition_spectral(adj(component_members,component_members));
    e = sparse(size(adj,1),1); e(component_members(find(idx_tmp==0)))=1; 
    components = [components, e]; isleaf = [isleaf,0]; isactive = [isactive, 1];
    e = sparse(size(adj,1),1); e(component_members(find(idx_tmp==1)))=1; 
    components = [components, e]; isleaf = [isleaf,0]; isactive = [isactive, 1];
    isactive(component_ind)=0;
end

idx = zeros(size(adj,1),1);
components = components(:,isleaf==1);
for i=1:size(components,2)
    idx(find(components(:,i)==1))=i;
end
return





% % this recursive structure is a bit slow
% function [idx, all_used_idx] = cut_mst(adj, data, current_idx, all_used_idx)
% 
% % these two variables are used to carry information during recursion
% if ~exist('current_idx')
%     current_idx = ones(size(adj,1));
% end
% if ~exist('all_used_idx')
%     all_used_idx=1;
% end
% threshold = 0.6;
% idx = zeros(size(adj,1),1);
% % check if the input mst satisfies stop criterion
% if avg_center_gene_corr(data)>threshold
%     idx = ones(size(adj,1),1)*current_idx;
%     return
% end
% % if the stop criterion is not met, we cut the mst into two parts
% [idx] = partition_longest_edge(adj);
% ind1 = find(idx==0); ind2 = find(idx==1);
% idx(ind1) = current_idx(1);
% idx(ind2) = max(all_used_idx)+1;
% all_used_idx = [all_used_idx, idx(ind2(1))];
% [idx(ind1), all_used_idx] = cut_mst(adj(ind1,ind1), data(ind1,:), idx(ind1), all_used_idx);
% [idx(ind2), all_used_idx] = cut_mst(adj(ind2,ind2), data(ind2,:), idx(ind2), all_used_idx);
% 
% return



function [idx] = partition_spectral(adj)
adj_TOM = TOM(TOM(adj));
L = diag(sum(adj_TOM))-adj_TOM; % identify elements in one piece
[V,D] = eigs(L,2,'sm');
V = V(:,1);
idx = double(V>=0);
if sum(idx)==0 || sum(idx)==size(adj,1)
   [idx] = partition_longest_edge(adj);
end
return



function [idx] = partition_longest_edge(adj)
% elimilate longest edge
[i,j] = find_matrix_top_element(adj); adj(i(1),j(1))=0; adj(j(1),i(1))=0;
% find pairs of connected nodes
[ind1,ind2] = find(triu(adj,1)~=0);
pairs = [ind1,ind2];
% identify elements in one piece
seed = i(1);
count =1;
while count<=length(seed)
    [a,b] = find(pairs==seed(count));
    seed = [seed; setdiff(unique(reshape(pairs(a,:),2*length(a),1)),seed)];
    count = count + 1;
end
idx = zeros(size(adj,1),1);
idx(seed)=1;
return
% [i,j] = find_matrix_top_element(adj);
% adj(i(1),j(1))=0; adj(j(1),i(1))=0;
% adj = double(adj~=0);
% adj = adj - diag(diag(adj)) + eye(size(adj));
% e = zeros(size(adj,1),1);
% e(i(1))=1; e(j(1))=-1;
% while sum(e==0)~=0
%     e = double(sign(adj*e));
% end
% idx = (e+1)/2;  % the returned value is either 0 or 1
return



