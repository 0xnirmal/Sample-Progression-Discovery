function [adj, idx] = get_adj_from_kmeans_phase1(data, rep, num_clusters)

idx = zeros(size(data,1),rep);
fprintf('Step 1: running multiple kmeans ... '); fprintf('%3d%%',0);
for i=1:size(idx,2)
    idx(:,i) = kmeans_phase1(data,num_clusters);
    fprintf('\b\b\b\b%3d%%',round(i/size(idx,2)*100)); drawnow;
end
fprintf('\b\b\b\bDone\n');

adj = zeros(size(data,1),size(data,1));
fprintf('Step 2: constructing adjacency ... '); fprintf('%3d%%',0);
for i=1:size(adj,1)-1
    for j=i+1:size(adj,1)
        adj(i,j) = sum(idx(i,:)==idx(j,:));
        adj(j,i) = adj(i,j);
    end
    adj(i,i) = max(adj(i,:));
    fprintf('\b\b\b\b%3d%%',round(i/size(adj,2)*100));drawnow
end
fprintf('\b\b\b\bDone\n');
