function idx = agglomerative_clustering(data, module_coherence_thres, merging_thres)
% each row if data matrix is one node (gene) to be clustered
% module_coherence_thres: target coherence threshold of a module, ~0.7, stopping criterion
% merging_thres: last heuristic step, merge similar modules, default is 0.9
%                if you want to skip this step, set merging_thres = 1;
% use the flow cell clustering algorithm, stop criterion is the coherence still use the final merge heuristic
% 
% data is per-gene-normalized to 0-mean 1-var

if ~exist('merging_thres')
    merging_thres = 0.9;
end

warning off;
data = per_gene_normalization(data); data = data./norm(data(1,:));

idx_genes = (1:size(data,1))';
module_means = data;
counter = 1;
while 1
    fold = max(2,round(max(idx_genes)/5000)); 
    fprintf('merging %d groups into about %d groups, round %d\n', max(idx_genes), ceil(max(idx_genes)/fold), counter);
    current_num_clusters = max(idx_genes);
    [module_means, idx_genes] = merge_genes_groups(module_means, idx_genes, fold, data,  module_coherence_thres);
    new_num_clusters= max(idx_genes);
    if current_num_clusters == new_num_clusters;
        break;
    end
    counter = counter + 1;
end 
idx = idx_genes;
module_coherence = zeros(max(idx),1);
for i=1:max(idx)
    module_coherence(i) = mean(module_means(i,:)*data(idx==i,:)');
end

fprintf('merging small pieces if they generate good avg_center_gene_corr ... \n')
fprintf('number of remaining_modules = %4d',max(idx));
% merge pieces if they generate good avg_center_gene_corr
c = module_means*module_means'; c = c - eye(size(c));
while max(max(c))>module_coherence_thres
    fprintf('\b\b\b\b%4d',max(idx));
    is_anything_merged=0;
    [pairs] = find_matrix_big_element(triu(c),module_coherence_thres);
    if isempty(pairs), break; end
    for i=1:size(pairs,1)
        data_tmp = data(idx==pairs(i,1) | idx==pairs(i,2),:);
        module_mean_tmp = mean(data_tmp,1); module_mean_tmp = module_mean_tmp/norm(module_mean_tmp);
        m = mean(module_mean_tmp*data_tmp');
        if m>=module_coherence_thres || c(pairs(i,1),pairs(i,2))>=merging_thres
            module_means(min(pairs(i,:)),:) = module_mean_tmp;
            c(min(pairs(i,:)),:) = module_means(min(pairs(i,:)),:)*module_means'; c(:,min(pairs(i,:))) = c(min(pairs(i,:)),:)'; c(min(pairs(i,:)),min(pairs(i,:)))=0;
            module_coherence(min(pairs(i,:)),:)=m;
            idx(idx==max(pairs(i,:))) = min(pairs(i,:));
            idx(idx>max(pairs(i,:))) = idx(idx>max(pairs(i,:)))-1;
            module_means(max(pairs(i,:)),:) = [];
            module_coherence(max(pairs(i,:)),:)=[];
            c(max(pairs(i,:)),:)=[];
            c(:,max(pairs(i,:)))=[];
            is_anything_merged=1;
            break;
        end
        c(pairs(i,1),pairs(i,2))=0; c(pairs(i,2),pairs(i,1))=0;
    end
    if is_anything_merged==0
        break;
    end
    drawnow
end
fprintf('\n');
% the above while loop
% (1) compute the corr between the module_means of all pairs of clusters (2)
% line 44 if the corr of a pair of module_means from two clusters is bigger than
% module_coherence_thres, it is worth looking into. For all the pairs that 
% are worth looking into, we order them according to the corr (3)
% according to the ordering, we look at these pairs one by one, see whether
% a pair can be merged, using two criteria: corr between > 0.9 or
% average coherence after merging > module_coherence_thres. (4) if merged,
% update c matrix (5) if not merged, this pairs' c is assigned to be zero,
% so that in the next round, we don't need to look at this pair again

return




function [module_means, idx] = merge_genes_groups(module_means, idx, fold, data,  module_coherence_thres)
% the genes grouped into give N groups by idx, 
% reduce the number of groups by half
% look at the groups one by one, 
% for each group, use single linkage to find the closet other group
% merge them together, and this other group is out of the game

isactive_module = ones(1,max(idx));
[module_size,X] = hist(idx,1:max(idx));
fprintf('merging in progress ... groups left: %7d  %7d',max(idx), max(idx));
iter = 1; total_num_modules = max(idx);
while sum(isactive_module)>1
%     module_ind = randsample(find(isactive_module==1),1);
    tmp = module_size; tmp(module_size==0)=Inf;  tmp(isactive_module==0) = Inf; 
    [dummy, module_ind] = min(tmp);
    module_center = module_means(module_ind,:);
    dist = module_center*(module_means');
    dist(module_ind) = 0;    
    [Y,I] = sort(dist,'descend'); 
    module_to_be_deleted = I(1:min(fold-1,end)); % this is the nearest PCA-cc-avglink modules
    if ~isempty(find(isactive_module(module_to_be_deleted)==0,1))
        module_to_be_deleted = setdiff(module_to_be_deleted, module_to_be_deleted(find(isactive_module(module_to_be_deleted)==0,1):end));
    end
    if isempty(module_to_be_deleted)
        isactive_module(module_ind)=0;   % if the nearest PCA-cc-avglink module is already examined, do not merge, continue on the while loop
    else
        for i=length(module_to_be_deleted):-1:1
            genes_ind = (idx==module_ind);
            for k=1:i
                genes_ind = (genes_ind | idx==module_to_be_deleted(k)); 
            end
            new_module_means = mean(data(genes_ind,:),1); new_module_means = new_module_means/norm(new_module_means);
            avg_coherence = abs(mean(new_module_means*(data(genes_ind,:)')));  % data is already normalized to 0-mean and norm 1 for each gene
            if avg_coherence<module_coherence_thres
                continue;
            else % merge them
                isactive_module(module_ind)=0;
                isactive_module(module_to_be_deleted(1:i))=0;
                idx(genes_ind) = module_ind;
                module_means(module_ind,:) = new_module_means;
                module_means(module_to_be_deleted(1:i),:) = 0;
                module_size(module_to_be_deleted(1:i))=0; module_size(module_ind) = sum(idx==module_ind);
                total_num_modules = total_num_modules - length(module_to_be_deleted((1:i)));
                break;
            end
        end
        isactive_module(module_ind)=0;
    end
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d  %7d',total_num_modules, sum(isactive_module)); drawnow;
%     iter = iter + 1; % iter is only serving for the drawnow in the following lines
%     if round(iter/10)==iter/10
%         drawnow;
%     end
end
fprintf('\n')
module_means = module_means(unique(idx),:);
idx = standardize_idx(idx);
return    



function idx_new = standardize_idx(idx)
% keep 0 as 0
% for positive integers, standardize to 1 2 3 4 ...
% if there are negative idx, the above will not work well

idx = idx(:); % make sure idx is a column vector
idx_new = zeros(size(idx));
[Y,I] = sort(idx,'ascend');
beginning_indicator_of_segments = [1;Y(1:end-1)~=Y(2:end)];
idx_new(I) = cumsum(beginning_indicator_of_segments);
if min(idx)==0
    idx_new = idx_new-1;
end




function [pairs] = find_matrix_big_element(matrix, threshold)
ind = find(matrix>=threshold);
j = ceil(ind/size(matrix,1));
i = ind - (j-1)*size(matrix,1);
pairs = [i,j];
[dummy,I] = sort(matrix(ind), 'descend');
pairs = pairs(I,:);
return


function data = per_gene_normalization(data)
% data = per_gene_normalization(data)

if size(data,2)<=1
    return
end

data = data - repmat(mean(data,2),1,size(data,2));
data = data./repmat(std(data')',1,size(data,2));
data(isnan(data))=0;









% function idx = agglomerative_clustering(data, module_coherence_thres, merging_thres)
% % each row if data matrix is one node (gene) to be clustered
% % module_coherence_thres: target coherence threshold of a module, ~0.7, stopping criterion
% % merging_thres: last heuristic step, merge similar modules, default is 0.9
% %                if you want to skip this step, set merging_thres = 1;
% % use the flow cell clustering algorithm, stop criterion is the coherence
% % still use the final merge heuristic
% 
% if ~exist('merging_thres')
%     merging_thres = 0.9;
% end
% 
% warning off;
% 
% idx_genes = 1:size(data,1);
% counter = 1;
% while 1
%     fold = max(2,round(max(idx_genes)/5000)); 
%     fprintf('merging %d groups into about %d groups, round %d\n', max(idx_genes), ceil(max(idx_genes)/fold), counter);
%     current_num_clusters = max(idx_genes);
%     idx_genes = merge_genes_to_groups(data, fold,idx_genes, module_coherence_thres);
%     new_num_clusters= max(idx_genes);
%     if current_num_clusters == new_num_clusters;
%         break;
%     end
%     counter = counter + 1;
% end 
% idx = idx_genes;
% 
% 
% 
% fprintf('merging small pieces if they generate good avg_center_gene_corr ... \n')
% fprintf('number of remaining_modules = %4d',max(idx));
% % merge pieces if they generate good avg_center_gene_corr
% [module_mean] = get_module_mean(data, idx);
% module_mean = per_gene_normalization(module_mean); module_mean = module_mean./norm(module_mean(1,:));
% c = module_mean*module_mean'; c = c - eye(size(c))*2;
% while max(max(c))>module_coherence_thres
%     fprintf('\b\b\b\b%4d',max(idx));
%     is_anything_merged=0;
%     [pairs] = find_matrix_big_element(triu(c),module_coherence_thres);
%     if isempty(pairs), break; end
%     for i=1:size(pairs,1)
%         m = avg_center_gene_corr(data(idx==pairs(i,1) | idx==pairs(i,2),:));
%         if m>=module_coherence_thres || c(pairs(i,1),pairs(i,2))>=merging_thres
%             % merge
%             module_mean(min(pairs(i,:)),:) = mean(data(idx==pairs(i,1) | idx==pairs(i,2),:),1);
%             module_mean(min(pairs(i,:)),:) = module_mean(min(pairs(i,:)),:)/norm(module_mean(min(pairs(i,:)),:));
%             c(min(pairs(i,:)),:) = module_mean(min(pairs(i,:)),:)*module_mean'; c(:,min(pairs(i,:))) = c(min(pairs(i,:)),:)'; c(min(pairs(i,:)),min(pairs(i,:)))=0;
%             idx(idx==max(pairs(i,:))) = min(pairs(i,:));
%             idx(idx>max(pairs(i,:))) = idx(idx>max(pairs(i,:))) - 1;
%             module_mean(max(pairs(i,:)),:)=[];
%             c(max(pairs(i,:)),:)=[];
%             c(:,max(pairs(i,:)))=[];
%             is_anything_merged=1;
%             break;
%         end
%         c(pairs(i,1),pairs(i,2))=0; c(pairs(i,2),pairs(i,1))=0;
%     end
%     if is_anything_merged==0
%         break;
%     end
%     drawnow
% end
% fprintf('\n');
% 
% return
% 
% 
% 
% function idx_new = merge_genes_to_groups(data, fold, idx, module_coherence_thres)
% % the genes grouped into give N groups by idx, 
% % reduce the number of groups by a factor of fold (min of fold is 2)
% % look at the groups one by one, 
% % for each group, use single linkage to find the closet other group
% % merge them together, and this other group is out of the game
% 
% isactive_module = ones(1,max(idx));
% [module_size,X] = hist(idx,1:max(idx)); 
% isactive_genes = ones(1,size(data,1));
% fprintf('merging in progress ... groups left: %7d  %7d',max(idx), max(idx));
% iter = 1; total_num_modules = max(idx);
% while sum(isactive_module)>1
% %     module_ind = randsample(find(isactive_module==1),1);
%     tmp = module_size; tmp(module_size==0)=Inf;  tmp(isactive_module==0) = Inf;
%     [dummy, module_ind] = min(tmp);
%     genes_in_module = find(idx==module_ind);
%     module_center = mean(data(genes_in_module,:),1);
%     dist = sum(abs(repmat(module_center,size(data,1),1) - data),2);
%     dist(idx==module_ind) = Inf;    % dist to everyone in my own group
%     [Y,I] = sort(dist,'ascend'); 
%     module_to_be_deleted = idx(I(1:fold-1));
%     first_inactive_module_ind = find(isactive_module(module_to_be_deleted)==0,1);
%     if ~isempty(first_inactive_module_ind) && first_inactive_module_ind==1, module_to_be_deleted=[]; end
%     if ~isempty(first_inactive_module_ind) && first_inactive_module_ind~=1, module_to_be_deleted=module_to_be_deleted(1:first_inactive_module_ind-1); end
%     if isempty(module_to_be_deleted)
%         isactive_module(module_ind)=0;
%         
%     else
%         
%         genes_to_be_merged_together = zeros(1,size(data,1)); genes_to_be_merged_together(idx == module_ind)=1;
%         for k=1:length(module_to_be_deleted), genes_to_be_merged_together(idx==module_to_be_deleted(k))=1; end
%         if avg_center_gene_corr(data(genes_to_be_merged_together==1,:))<module_coherence_thres
%             isactive_module(module_ind)=0;
%         else
%             isactive_module(module_to_be_deleted)=0;
%             isactive_module(module_ind)=0;
%             for i=1:length(module_to_be_deleted)
%                 idx(idx==module_to_be_deleted(i)) = module_ind;  % merge
%             end
%             module_size(module_to_be_deleted)=0; module_size(module_ind) = sum(idx==module_ind);
%             isactive_sample(idx==module_ind)=0;
%             total_num_modules = total_num_modules - length(module_to_be_deleted);
%         end
%         
%     end    
%     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d  %7d',total_num_modules, sum(isactive_module)); drawnow;
%     drawnow;
% end
% fprintf('\n')
% idx_new = standardize_idx(idx);
% return    
% 
