function [idx,PCAs,module_coherence] = agglomerative_clustering_absCC_PCAlink(data, module_coherence_thres, merging_thres)
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
PCAs = data;
counter = 1;
while 1
    fold = max(2,round(max(idx_genes)/5000)); 
    fprintf('merging %d groups into about %d groups, round %d\n', max(idx_genes), ceil(max(idx_genes)/fold), counter);
    current_num_clusters = max(idx_genes);
    [PCAs, idx_genes] = merge_genes_groups(PCAs, idx_genes, fold, data,  module_coherence_thres);
    new_num_clusters= max(idx_genes);
    if current_num_clusters == new_num_clusters;
        break;
    end
    counter = counter + 1;
end 
idx = idx_genes;
module_coherence = zeros(max(idx),1);
for i=1:max(idx)
    module_coherence(i) = mean(abs(PCAs(i,:)*data(idx==i,:)'));
end

fprintf('merging small pieces if they generate good avg_center_gene_corr ... \n')
fprintf('number of remaining_modules = %4d',max(idx));
% merge pieces if they generate good avg_center_gene_corr
c = abs(PCAs*PCAs'); c = c - eye(size(c));
while max(max(c))>module_coherence_thres
    fprintf('\b\b\b\b%4d',max(idx));
    is_anything_merged=0;
    [pairs] = find_matrix_big_element(triu(c),module_coherence_thres);
    if isempty(pairs), break; end
    for i=1:size(pairs,1)
        data_tmp = data(idx==pairs(i,1) | idx==pairs(i,2),:);
        [U,S,V] = svd(data_tmp); PCA_tmp = V(:,1)';
        m = mean(abs(PCA_tmp*data_tmp'));
        if m>=module_coherence_thres || c(pairs(i,1),pairs(i,2))>=merging_thres
            PCAs(min(pairs(i,:)),:) = PCA_tmp;
            c(min(pairs(i,:)),:) = abs(PCAs(min(pairs(i,:)),:)*PCAs'); c(:,min(pairs(i,:))) = c(min(pairs(i,:)),:)'; c(min(pairs(i,:)),min(pairs(i,:)))=0;
            module_coherence(min(pairs(i,:)),:)=m;
            idx(idx==max(pairs(i,:))) = min(pairs(i,:));
            idx(idx>max(pairs(i,:))) = idx(idx>max(pairs(i,:)))-1;
            PCAs(max(pairs(i,:)),:) = [];
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
% (1) compute the abs corr between the PCAs of all pairs of clusters (2)
% line 44 if the abs corr of a pair of PCAs from two clusters is bigger than
% module_coherence_thres, it is worth looking into. For all the pairs that 
% are worth looking into, we order them according to the abs corr (3)
% according to the ordering, we look at these pairs one by one, see whether
% a pair can be merged, using two criteria: abs corr between > 0.9 or
% average coherence after merging > module_coherence_thres. (4) if merged,
% update c matrix (5) if not merged, this pairs' c is assigned to be zero,
% so that in the next round, we don't need to look at this pair anymore

return




function [PCAs, idx] = merge_genes_groups(PCAs, idx, fold, data,  module_coherence_thres)
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
    module_center = PCAs(module_ind,:);
    dist = abs(module_center*(PCAs'));
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
            [U,S,V] = svd(data(genes_ind,:));
            newPCA = V(:,1)';
            avg_coherence = mean(abs(newPCA*(data(genes_ind,:)')));  % data is already normalized to 0-mean and norm 1 for each gene
            if avg_coherence<module_coherence_thres
                continue;
            else % merge them
                isactive_module(module_ind)=0;
                isactive_module(module_to_be_deleted(1:i))=0;
                idx(genes_ind) = module_ind;
                PCAs(module_ind,:) = newPCA;
                PCAs(module_to_be_deleted(1:i),:) = 0;
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
PCAs = PCAs(unique(idx),:);
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