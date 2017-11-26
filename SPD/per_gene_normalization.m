function data = per_gene_normalization(data)
% data = per_gene_normalization(data)

if size(data,2)<=1
    return
end

data = data - repmat(mean(data,2),1,size(data,2));
data = data./repmat(std(data')',1,size(data,2));
data(isnan(data))=0;
