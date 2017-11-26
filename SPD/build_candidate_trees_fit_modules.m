function [candidate_trees, p_value] = build_candidate_trees_fit_modules(data,idx, permutation_iter)
% this is\was part of the file analyze_progression
% when making the GUI, it seems good to break that file into parts
% later, I may use this function to replace the corresponding part in 
% analyze_progression, to make it look cleaner


iter = permutation_iter;
candidate_trees = zeros(size(data,2),size(data,2),max(idx));
for i=1:max(idx), 
    module_size(i) = sum(idx==i); 
    [adj1,adj2,cost] = mst(data(idx==i,:)');
    candidate_trees(:,:,i) = adj1;
    drawnow;
end

fprintf('\n\n finished building candidate trees for %d modules \n\n', max(idx));

% % for each tree and each module, pair, find out how they fit each other
p_value = [];
fprintf('fitting candidate trees and modules, pairs %4d, %4d', 0,0);
for j=1:max(idx)
    distance_matrix = get_row_distance(data(idx==j,:)');
    for i=1:size(candidate_trees,3)
        progression_tree = candidate_trees(:,:,i);
        fprintf('\b\b\b\b\b\b\b\b\b\b%4d, %4d', i,j); drawnow;
        score = fit_data_mst(progression_tree, data(idx==j,:)); % fit_data_mst(progression_tree, data(idx==j,:), [],distance_matrix)

        null_scores = zeros(1,iter);
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
