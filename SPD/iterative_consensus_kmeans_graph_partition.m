function idx = iterative_consensus_kmeans(data, kmeans_iter, module_coherence_thres, merging_thres)
% each row if data matrix is one node (gene) to be clustered
% kmeans_iter: num of kmenas iterations in the consensus algorithm
% module_coherence_thres: target coherence threshold of a module, ~0.7, stopping criterion
% merging_thres: last heuristic step, merge similar modules, default is 0.9
%                if you want to skip this step, set merging_thres = 1;


if ~exist('merging_thres')
    merging_thres = 0.9;
end

warning off;

components = sparse(ones(size(data,1),1));
if avg_center_gene_corr(data(components(:,1)==1,:))>module_coherence_thres
    isleaf = 1; 
    isactive=0;
else
    isleaf = 0;
    isactive=1;
end

fprintf('performing top-down consensus kmeans ... \n')
fprintf('number of processed intermediate nodes / total intermediate nodes %4d /%4d',0,1);
idx_tmp = zeros(size(data,1),1);
while sum(isactive==1)~=0
    component_ind = find(isactive==1,1);
    fprintf('\b\b\b\b\b\b\b\b\b\b%4d /%4d',component_ind-1, length(isactive));
    component_members = find(components(:,component_ind)==1);
    if length(component_members)==1
        isleaf(component_ind) = 1; 
        isactive(component_ind)=0;
        continue;
    end
    if avg_center_gene_corr(data(component_members,:))>module_coherence_thres
        isleaf(component_ind) = 1; 
        isactive(component_ind)=0;
        continue;
    end
    [idx_tmp] = consensus_2means_partition(data(component_members,:), kmeans_iter);
    e = sparse(size(data,1),1); e(component_members(find(idx_tmp==0)))=1; 
    components = [components, e]; isleaf = [isleaf,0]; isactive = [isactive, 1];
    e = sparse(size(data,1),1); e(component_members(find(idx_tmp==1)))=1; 
    components = [components, e]; isleaf = [isleaf,0]; isactive = [isactive, 1];
    if sum(components(:,end-1))==0 || sum(components(:,end))==0
        1;
    end
    isactive(component_ind)=0;
    drawnow
end
components_all = components;
fprintf('\nDone\n');

idx = zeros(size(data,1),1);
components = components(:,isleaf==1);
for i=1:size(components,2)
    idx(find(components(:,i)==1))=i;
end



fprintf('merging small pieces if they generate good avg_center_gene_corr ... \n')
fprintf('number of remaining_modules = %4d',max(idx));
% merge pieces if they generate good avg_center_gene_corr
[module_mean] = get_module_mean(data, idx);
module_mean = per_gene_normalization(module_mean); module_mean = module_mean./norm(module_mean(1,:));
c = module_mean*module_mean'; c = c - eye(size(c))*2;
while max(max(c))>module_coherence_thres
    fprintf('\b\b\b\b%4d',max(idx));
    is_anything_merged=0;
    [pairs] = find_matrix_big_element(triu(c),module_coherence_thres);
    if isempty(pairs), break; end
    for i=1:size(pairs,1)
        m = avg_center_gene_corr(data(idx==pairs(i,1) | idx==pairs(i,2),:));
        if m>=module_coherence_thres || c(pairs(i,1),pairs(i,2))>=merging_thres
            % merge
            module_mean(min(pairs(i,:)),:) = mean(data(idx==pairs(i,1) | idx==pairs(i,2),:),1);
            module_mean(min(pairs(i,:)),:) = module_mean(min(pairs(i,:)),:)/norm(module_mean(min(pairs(i,:)),:));
            c(min(pairs(i,:)),:) = module_mean(min(pairs(i,:)),:)*module_mean'; c(:,min(pairs(i,:))) = c(min(pairs(i,:)),:)'; c(min(pairs(i,:)),min(pairs(i,:)))=0;
            idx(idx==max(pairs(i,:))) = min(pairs(i,:));
            idx(idx>max(pairs(i,:))) = idx(idx>max(pairs(i,:))) - 1;
            module_mean(max(pairs(i,:)),:)=[];
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

% 
% fprintf('merging small pieces if they generate good avg_center_gene_corr ... \n')
% fprintf('number of remaining_modules = %4d',max(idx));
% % merge pieces if they generate good avg_center_gene_corr
% [module_mean] = get_module_mean(data, idx);
% module_mean = per_gene_normalization(module_mean); module_mean = module_mean./norm(module_mean(1,:));
% c = module_mean*module_mean'; c = c - eye(size(c))*2;
% while max(max(c))>module_coherence_thres
%     fprintf('\b\b\b\b%4d',max(idx));
%     [pairs] = find_matrix_big_element(triu(c),module_coherence_thres);
%     if isempty(pairs), break; end
%     tmp_merge_corr = zeros(size(pairs,1),1);
%     for i=1:size(pairs,1)
%         tmp_merge_corr(i) = avg_center_gene_corr(data(idx==pairs(i,1) | idx==pairs(i,2),:));
%     end
%     [m,mind] = max(tmp_merge_corr);
%     if m>=module_coherence_thres || c(pairs(mind,1),pairs(mind,2))>=merging_thres
%         idx(idx==pairs(mind,2)) = pairs(mind,1);
%         idx = standardize_idx(idx);
%         [module_mean] = get_module_mean(data, idx);
%         module_mean = per_gene_normalization(module_mean); module_mean = module_mean./norm(module_mean(1,:));
%         c = module_mean*module_mean'; c = c - eye(size(c))*2;
% %         module_mean = normalize(module_mean')';
%     else
%         break;
%     end
%     drawnow
% end
% fprintf('\n');

return



function [idx] = consensus_2means_partition(data, kmeans_iter)
% the elements of returned idx are either 0 or 1

if size(data,1)==2
    idx = [0;1]; return
end

num_clusters=2;
idx_iter = zeros(size(data,1),kmeans_iter);
fprintf(' ... %3d%%',0); drawnow;
for i=1:kmeans_iter
try
    idx_iter(:,i) = kmeans_phase1(data,num_clusters);
catch
    1
end
    fprintf('\b\b\b\b%3d%%',round(i/kmeans_iter*100)); drawnow;
%     drawnow;
end
fprintf('\b\b\b\b\b\b\b\b\b'); drawnow;
% idx contains the kmeans clustering results,
% if two genes are similar, then the two rows are similar, 
% therefore, we can simply treat one row of idx as the expression of a gene
% and perform another round of clustering
d = pdist(idx_iter,'hamming');
[m, mind] = max(d);
for i=1:size(data,1)-1, if mind<=size(data,1)-i, break; else, mind = mind - (size(data,1)-i); end; end
% the above two lines finds the two genes that are most different from each
% other in the multiple kmeans results: idx_iter matrix
idx = kmeans_phase1(idx_iter, 2, 'Start', idx_iter([i,i+mind],:))-1;
return



