

function [avg_corr_value] = avg_center_gene_corr(data)

data = per_gene_normalization(data); data = data./norm(data(1,:));
m = mean(data,1);
m = m./norm(m);
avg_corr_value = mean(m*data');

% % average gene to center correlation
% m = mean(data,1);
% 
% % remove the mean component
% m = m-mean(m);
% data = data - repmat(mean(data,2),1,size(data,2));
% 
% data_norm = sqrt(sum(data.^2,2)); data_norm(data_norm==0) = 1; % this handles the case if one row of data (one gene) have expression equals to 0 across all samples
% m_norm = sqrt(sum(m.^2)); m_norm(m_norm==0) = 1;
% corr_values = data*m'./data_norm/m_norm;
% avg_corr_value = mean(corr_values);
