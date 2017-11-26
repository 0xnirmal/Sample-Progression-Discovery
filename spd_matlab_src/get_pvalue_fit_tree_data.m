function [p_value, score, null_scores] = get_pvalue_fit_tree_data(progression_tree,data)

distance_matrix = get_row_distance(data'); 
score = fit_data_mst(progression_tree, data); % fit_data_mst(progression_tree, data(idx==j,:), [],distance_matrix)

iter = 1000; null_scores = zeros(1,iter);
% fprintf('performing permutation to fit tree and data ... %5d',0);
for k=1:iter, 
    null_perm = randsample(1:size(data,2),size(data,2));
    null_scores(k) = sum(sum(triu(progression_tree).*distance_matrix(null_perm,null_perm)))/sum(sum(triu(progression_tree))); 
%     fprintf('\b\b\b\b\b%5d',k);
    % fitting_score_null(k) = fit_data_mst(progression_tree, data(idx==j,null_perm), [],distance_matrix(null_perm,null_perm)); 
end
% fprintf('\n');
p_value =  sum(score>=null_scores)/iter;  % normcdf(score, mean(fitting_score_null), std(fitting_score_null));