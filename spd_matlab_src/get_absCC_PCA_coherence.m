function m = get_absCC_PCA_coherence(data, idx)
% m = get_absCC_PCA_coherence(data, idx)
% (idx is the clustering result of genes/rows)

data = per_gene_normalization(data); data = data./norm(data(1,:));

m = zeros(1,max(idx));
for i=1:max(idx)
    data_tmp = data(idx==i,:);
    [U,S,V] = svd(data_tmp); PCA_tmp = V(:,1)';
    m(i) = mean(abs(PCA_tmp*data_tmp'));
end

return
