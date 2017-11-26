

function [avg_corr_value] = avg_center_gene_corr(data)
% average gene to center correlation
m = mean(data,1);
data_norm = sqrt(sum(data.^2,2)); data_norm(data_norm==0) = 1; % this handles the case if one row of data (one gene) have expression equals to 0 across all samples
m_norm = sqrt(sum(m.^2)); m_norm(m_norm==0) = 1;
corr_values = data*m'./data_norm/m_norm;
avg_corr_value = mean(corr_values);
