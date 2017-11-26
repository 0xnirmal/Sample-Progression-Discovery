function display_mst(tree_for_embedding)

tree_backup = tree_for_embedding;
tree_for_embedding(triu(tree_for_embedding,1)==1) = 1 + rand(sum(sum(triu(tree_for_embedding,1)==1)),1)*10;
tree_for_embedding = (tree_for_embedding + tree_for_embedding')/2;
[positions] = hd_embedding(tree_for_embedding,5);
hold off; [coeff] = multiclass_pcaplot(positions); hold on;
pairs = find_matrix_big_element(triu(tree_backup,1),1);
for k=1:size(pairs,1), line(coeff(1,pairs(k,:)),coeff(2,pairs(k,:)),'color','g'); end
plot(coeff(1,:),coeff(2,:),'o','markersize',3,'markerfacecolor','b');
for k=1:size(coeff,2), text(coeff(1,k)+1,coeff(2,k),num2str(k)); end  
