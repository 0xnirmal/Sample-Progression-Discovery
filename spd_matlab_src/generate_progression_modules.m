function progression_modules = generate_progression_modules(data,idx, merging_pvalue, permutation_iter)
% this is\was part of the file analyze_progression
% when making the GUI, it seems good to break that file into parts
% later, I may use this function to replace the corresponding part in 
% analyze_progression, to make it look cleaner

p_threshold = merging_pvalue;
iter = permutation_iter;

% initialize progression_modules, each progression module contains one gene
% module, the progression modules are sorted according to their size (how
% many genes are there)
gene_module_mean = get_module_mean(data,idx);
for i=1:max(idx)
    progression_modules(i).gene_modules = i; 
    progression_modules(i).gene_modules_mean = gene_module_mean(i,:);
    progression_modules(i).gene_expr = data(idx==i,:); %./sum(idx==i);
    progression_module_size(i) = size(progression_modules(i).gene_expr,1);
    [adj1,adj2,cost] = mst(progression_modules(i).gene_expr');
    progression_modules(i).adj1 = adj1;
    progression_modules(i).adj2 = adj2;
    progression_modules(i).node_position = [];
end
[Y,I] = sort(progression_module_size,'descend');
progression_modules = progression_modules(I);  % sorting progression modules according to number of genes in them
isactive = ones(size(progression_modules));
return
% use the top progression_module (every gene in it) to build a tree, see
% whether other progression_modules fit well with this tree, if there are
% any, merge the best into the top progression module, rebuild the tree, to
% it again. when no more merge is acceptable (there is a p-threshold for
% merging), put it aside as one progression_module to report, the second
% one becomes the top module, do it all over again
while sum(isactive)~=0
    ind = find(isactive==1,1);
    p_values = ones(size(progression_modules));
    fprintf('comparing top module %d vs all %d modules ... %4d',ind, length(p_values), 0);
    for i=1:length(progression_modules)
        if isactive(i)==0 || i==ind, continue; end
%         [p_values(i), score, null_scores] = get_pvalue_fit_tree_data(progression_modules(ind).adj1,progression_modules(i).gene_expr);
        progression_tree = progression_modules(ind).adj1; distance_matrix = get_row_distance(progression_modules(i).gene_expr');
        score = fit_data_mst(progression_tree, progression_modules(i).gene_expr);
        null_scores = zeros(1,iter);
        non_zero_ind = find(triu(progression_tree)~=0); denominator = sum(progression_tree(non_zero_ind));
        for k=1:iter, 
            [dummy, null_perm] = sort(rand(1,size(progression_modules(i).gene_expr,2)));
            distance_matrix_tmp = distance_matrix(null_perm,null_perm);
            null_scores(k) = sum(progression_tree(non_zero_ind).*distance_matrix_tmp(non_zero_ind))/denominator;
        end
        p_values(i) = sum(score>=null_scores)/iter; 

        fprintf('\b\b\b\b%4d',i); drawnow;
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