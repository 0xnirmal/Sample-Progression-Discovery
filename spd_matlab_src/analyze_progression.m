function [idx, progression_modules] = identify_progression(data,probe_names,filename,gene_cluster_t, module_cluster_t1, module_cluster_t2)

if ~exist('gene_cluster_t'), gene_cluster_t = 0.7; end
if ~exist('module_cluster_t1'), module_cluster_t1 = 0.8; end
if ~exist('module_cluster_t2'), module_cluster_t2 = 0.9; end


if ~exist('filename');
    filename = 'tmp_progression_results.mat';
end

% % concensus kmeans clustering of genes, results in "idx"
[adj] = get_adj_from_kmeans_phase1(data,200,100);
[perm_r,perm_c] = HCC_heatmap(adj);
save(filename,'data','probe_names','adj','perm_r');
% % load(filename)
[idx] = cut_corr_matrix(adj/200,data,gene_cluster_t,0.9);
figure(3); [Y,I] = sort(idx); imagesc(adj(I,I)); drawnow% this line shows that my code has done an OK job in clustering
save(filename,'idx','gene_cluster_t','-append');
% % ss=[];for i=1:max(idx), ss(i) = avg_center_gene_corr(data(idx==i,:)); end
fprintf('\n\n clustering of genes finishsed, resulting in %d modules of a total of %d genes\n\n', max(idx), size(data,1));

% % generate progression_modules
gene_module_mean = get_module_mean(data,idx);
for i=1:max(idx)
    progression_modules(i).gene_modules = i; 
    progression_modules(i).gene_modules_mean = gene_module_mean(i,:);
    progression_modules(i).gene_expr = data(idx==i,:); %./sum(idx==i);
    progression_module_size(i) = size(progression_modules(i).gene_expr,1);
    [adj1,adj2,cost] = mst(progression_modules(i).gene_expr');
    progression_modules(i).adj1 = adj1;
    progression_modules(i).adj2 = adj2;
end
[Y,I] = sort(progression_module_size,'descend');
progression_modules = progression_modules(I); isactive = ones(size(progression_modules));
p_threshold = 0.001;
while sum(isactive)~=0
    ind = find(isactive==1,1);
    p_values = ones(size(progression_modules));
    fprintf('comparing top module %d vs all %d modules ... %4d',ind, length(p_values), 0);
    for i=1:length(progression_modules)
        if isactive(i)==0 || i==ind, continue; end
%         [p_values(i), score, null_scores] = get_pvalue_fit_tree_data(progression_modules(ind).adj1,progression_modules(i).gene_expr);
        progression_tree = progression_modules(ind).adj1; distance_matrix = get_row_distance(progression_modules(i).gene_expr');
        score = fit_data_mst(progression_tree, progression_modules(i).gene_expr);
        iter = 1000; null_scores = zeros(1,iter);
        non_zero_ind = find(triu(progression_tree)~=0); denominator = sum(progression_tree(non_zero_ind));
        for k=1:iter, 
            [dummy, null_perm] = sort(rand(1,size(progression_modules(i).gene_expr,2)));
            distance_matrix_tmp = distance_matrix(null_perm,null_perm);
            null_scores(k) = sum(progression_tree(non_zero_ind).*distance_matrix_tmp(non_zero_ind))/denominator;
        end
        p_values(i) = sum(score>=null_scores)/iter; 

        fprintf('\b\b\b\b%4d',i);
    end
    fprintf('\n');
    [m,mind] = min(p_values);
    fprintf('best fit p-value %1.4f\n',m);
    if m>p_threshold,
        isactive(ind)=0; continue;
    end
    progression_modules(ind).gene_modules = [progression_modules(ind).gene_modules; progression_modules(mind).gene_modules]; 
    progression_modules(ind).gene_modules'
    progression_modules(ind).gene_modules_mean = [progression_modules(ind).gene_modules_mean; progression_modules(mind).gene_modules_mean];
    progression_modules(ind).gene_expr = [progression_modules(ind).gene_expr;progression_modules(mind).gene_expr];
    [adj1,adj2,cost] = mst(progression_modules(ind).gene_expr');
    progression_modules(ind).adj1 = adj1;
    progression_modules(ind).adj2 = adj2;
    progression_modules(mind)=[];
    isactive(mind)=[];
end
save(filename,'progression_modules','-append');



% 
% % for each module define in "idx", compute a candidate minimum spanning tree
candidate_trees = zeros(size(data,2),size(data,2),max(idx));
for i=1:max(idx), 
    module_size(i) = sum(idx==i); 
    [adj1,adj2,cost] = mst(data(idx==i,:)');
    candidate_trees(:,:,i) = adj1;
end
save(filename,'candidate_trees','module_size','-append');

fprintf('\n\n finished building candidate trees for %d modules \n\n', max(idx));


% % for each tree and each module, pair, find out how they fit each other
p_value = [];
fprintf('fitting candidate trees and modules, pairs %4d, %4d', 0,0);
for j=1:max(idx)
    distance_matrix = get_row_distance(data(idx==j,:)');
    for i=1:size(candidate_trees,3)
        progression_tree = candidate_trees(:,:,i);
        fprintf('\b\b\b\b\b\b\b\b\b\b%4d, %4d', i,j);
        score = fit_data_mst(progression_tree, data(idx==j,:)); % fit_data_mst(progression_tree, data(idx==j,:), [],distance_matrix)

        iter = 1000; null_scores = zeros(1,iter);
        non_zero_ind = find(triu(progression_tree)~=0); denominator = sum(progression_tree(non_zero_ind));
        for k=1:iter, 
            [dummy, null_perm] = sort(rand(1,size(data,2)));
            distance_matrix_tmp = distance_matrix(null_perm,null_perm);
            null_scores(k) = sum(progression_tree(non_zero_ind).*distance_matrix_tmp(non_zero_ind))/denominator;
        end
        p_value(i,j) =  sum(score>=null_scores)/iter;
        
%         iter = 1000; fitting_score_null = zeros(1,iter);
%         for k=1:iter, 
%             null_perm = randsample(1:size(data,2),size(data,2));
%             fitting_score_null(k) = sum(sum(triu(progression_tree).*distance_matrix(null_perm,null_perm)))/sum(sum(triu(progression_tree))); 
%             % fitting_score_null(k) = fit_data_mst(progression_tree, data(idx==j,null_perm), [],distance_matrix(null_perm,null_perm)); 
%         end
%         p_value(i,j) =  sum(score>=fitting_score_null)/iter;  % normcdf(score, mean(fitting_score_null), std(fitting_score_null));
    end
end
save(filename,'p_value','-append');


% % construct an adjacency of modules, if two modules both fit one/more
% % candidate trees, these two modules are more strongly connected
p_threshold = 0.002;
adj_modules = zeros(size(candidate_trees,3));
for i=1:size(p_value,1)
    ind_tmp = find(p_value(i,:)<=p_threshold);
    adj_modules(ind_tmp,ind_tmp)=adj_modules(ind_tmp,ind_tmp) + 1;
end
figure(4); perm_modules = HCC_heatmap(adj_modules); axis off
for i=1:length(perm_modules)
    text(i,-1,num2str(perm_modules(i)));
    text(i,length(perm_modules)+2,num2str(perm_modules(i)));
end

fprintf('\n\n')
