function [idx, components_all] = cut_corr_matrix(corr_matrix,data,threshold, threshold2)

if ~exist('threshold2')
    threshold2 = 0.9;
end

warning off;

adj = corr_matrix;

% components = sparse(ones(size(adj,1),1));
% if avg_center_gene_corr(data(components(:,1)==1,:))>threshold
%     isleaf = 1; 
%     isactive=0;
% else
%     isleaf = 0;
%     isactive=1;
% end
fprintf('initialization by identify disconnected components ... \n');
[components] = initialization_cut_corr_matrix(adj); 
isleaf = zeros(1,size(components,2)); isactive = zeros(1,size(components,2));
for i=1:size(components,2)
    if sum(components(:,i))==1 || avg_center_gene_corr(data(components(:,1)==1,:))>threshold
        isleaf(i) = 1; isactive(i) = 0;
    else
        isleaf(i) = 0; isactive(i) = 1;
    end
end

fprintf('performing top-down spectral divide ... \n')

idx_tmp = zeros(size(adj,1),1);
while sum(isactive==1)~=0
    component_ind = find(isactive==1,1);
%     [component_ind, length(isactive)]
    component_members = find(components(:,component_ind)==1);
    [subcomponents] = initialization_cut_corr_matrix(adj(component_members,component_members));
    if size(subcomponents,2)~=1  % which means that this component is already disconnected
        isleaf(component_ind) = 0; 
        isactive(component_ind)=0;
        tmp = zeros(size(components,1),size(subcomponents,2)); tmp(component_members,:) = subcomponents;
        components = [components, tmp];
        isleaf = [isleaf, zeros(1,size(subcomponents,2))];
        isactive = [isactive, ones(1,size(subcomponents,2))];
        continue;
    end
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
    if sum(components(:,end-1))==0 || sum(components(:,end))==0
        1;
    end
    isactive(component_ind)=0;
    drawnow
end
components_all = components;

idx = zeros(size(adj,1),1);
components = components(:,isleaf==1);
for i=1:size(components,2)
    idx(find(components(:,i)==1))=i;
end


fprintf('merging small pieces if they generate good avg_center_gene_corr ... \n')

% merge pieces if they generate good avg_center_gene_corr
[module_mean] = get_module_mean(data, idx);
module_mean = per_gene_normalization(module_mean); module_mean = module_mean./norm(module_mean(1,:));
c = module_mean*module_mean'; c = c - eye(size(c))*2;
while max(max(c))>threshold
%     max(idx)
    [pairs] = find_matrix_big_element(triu(c),threshold);
    if isempty(pairs), break; end
    tmp_merge_corr = zeros(size(pairs,1),1);
    for i=1:size(pairs,1)
        tmp_merge_corr(i) = avg_center_gene_corr(data(idx==pairs(i,1) | idx==pairs(i,2),:));
    end
    [m,mind] = max(tmp_merge_corr);
    if m>=threshold || c(pairs(mind,1),pairs(mind,2))>=threshold2
        idx(idx==pairs(mind,2)) = pairs(mind,1);
        idx = standardize_idx(idx);
        [module_mean] = get_module_mean(data, idx);
        module_mean = per_gene_normalization(module_mean); module_mean = module_mean./norm(module_mean(1,:));
        c = module_mean*module_mean'; c = c - eye(size(c))*2;
%         module_mean = normalize(module_mean')';
    else
        break;
    end
    drawnow
end


return




function [idx] = partition_spectral(adj)
L = sparse(diag(sum(adj)))-sparse(adj); % identify elements in one piece
if sum(sum(abs(L)))<1e-10 % L is virtually zero matrix
    idx = zeros(size(adj,1),1);
    idx(randsample(1:end,round(length(idx)/2)))=1;
else
    OPTS.disp = 0;
    [V,D] = eigs(L,2,'sm',OPTS);
    [dummy,D_ind] = max(abs(diag(D))); V = V(:,D_ind);
    if isnan(dummy)
        [V,D] = eig(full(L)); 
        [Y,I] = sort(abs(diag(D))); D_ind = I(2); V = V(:,D_ind);
    end
    idx = double(V>=0);
    if sum(idx)==0 || sum(idx)==size(adj,1)
       disp('error');
    end
end
return



function [components] = initialization_cut_corr_matrix(adj)
% first find disconnected components
adj_matrix = adj>0;
for i=1:size(adj_matrix,1), adj_matrix(i,i)=1; end
components=[];
is_assigned = zeros(size(adj_matrix,1),1);
is_assigned = full(is_assigned);
while sum(is_assigned==0)~=0
%     [sum(is_assigned),length(is_assigned)]
    e = zeros(size(adj_matrix,1),1);
    e(find(is_assigned==0,1))=1;
    while sum(e~=(adj_matrix*e>0))~=0
        e = double((adj_matrix*e)>0);
    end
    components = [components,e];
    is_assigned(find(e==1))=1;
end
return
