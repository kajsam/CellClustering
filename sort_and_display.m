function [idx, Zet] = sort_and_display(X, Zet, A, Pset, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 12th 2018

% Input:        X - binary matrix
%               Zsets - the sets
%               A - fisher matrix

[n,d] = size(X);
K = size(Zet,2);

[~,idx] = sort(sum(Zet),'descend');
Zet = Zet(:,idx);

figure(fig_nr), imagesc(Zet), colormap(gray), title('Log likelihood cell clusters')
figure, imagesc(A), colormap(gray), title('Log likelihood Fisher matrix')

id = zeros(1,K);
for k = 1: K
  id(k) = find(Zet(:,k),1);
end

% Number of possible gene classes and their names
name_gene_class_top = cell(1);
name_gene_class_bottom = cell(1);
t = 0; 
for k = 1: K
  c = combnk(1:K,k);
  for i = 1: size(c,1)
    t = t+1;
    name_gene_class_top{t} = num2str(c(i,:));
    name_gene_class_top{t}
    name_gene_class_bottom{t} = num2str(c(i,:));
  end
end
t = t+1;
name_gene_class_top{t} = '0';
name_gene_class_bottom{t} = '0';

for i = 1:2:t
  name_gene_class_top{i} = ' ';
end
for i = 2:2:t
  name_gene_class_bottom{i} = ' ';
end

class = t*ones(1,d);
for j = 1: d
  class_name = 0;
  for k = 1 : K
    c = combnk(1:K,k);
    for i = 1: size(c,1)
      class_name = class_name+1;
      if all(A(id(c(i,:)),j))  
        class(j) = class_name;
      end
    end
  end
end

[C,~,ic] = unique(class);
a_counts = accumarray(ic,1);
[C', a_counts]

idx_c = cell(1,t);
for c = 1: t
  idx_c{c} = find(class == c);
end

for c = 1: t
  [~,idx] = sort(sum(X(:,idx_c{c})),'descend');
  idx_c{c} = idx_c{c}(idx);
end

idx = cell2mat(idx_c);

figure(fig_nr+1), imagesc(X(:,idx)), colormap(gray)
xticks(cumsum([1; a_counts(1:end-1)]))
xticklabels(name_gene_class_bottom)
% ax1=gca;
% ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
% set(ax2, 'XAxisLocation', 'top');
% % set the same Ticks on ax2 as on ax1;
% set(ax2, 'XTick', get(ax1, 'XTick'));
% % % Set the x-tick labels for the second axes
% set(ax2, 'XTickLabel', name_gene_class_top);

idx_c = cell(1,t);
for c = 1: t
  idx_c{c} = find(class == c);
end

for c = 1: t
  [~,idxp] = sort(Pset(:,idx_c{c}),'ascend');
  idx_c{c} = idx_c{c}(idxp);
end

idxp = cell2mat(idx_c);

% figure(fig_nr+2), imagesc(X(:,idxp)), colormap(gray)
% xticks(cumsum(a_counts(1:end-1)))
% xticklabels(name_gene_class)



    