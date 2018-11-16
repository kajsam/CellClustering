function [idx, Zet, ticks, ticklabels] = sort_and_display(X, Zet, A, Pset, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 12th 2018

% Input:        X - binary matrix
%               Zsets - the sets
%               A - fisher matrix

[n,d] = size(X);
K = size(Zet,2);

% [~,idx] = sort(sum(Zet),'descend');
% Zet = Zet(:,idx);

figure(fig_nr), imagesc(Zet), colormap(gray), title('Log likelihood cell clusters')
figure(fig_nr+1), imagesc(A), colormap(gray), title('Log likelihood Fisher matrix')

id = zeros(1,K);
for k = 1: K
  id(k) = find(Zet(:,k),1);
end

% Number of possible gene classes and their names
name_gene_class = cell(1);
%name_gene_class_bottom = cell(1);
t = 0; 
for k = 1: K
  c = combnk(1:K,k);
  for i = 1: size(c,1)
    t = t+1;
    name_gene_class{t} = num2str(c(i,:));
    % name_gene_class_top{t}
    %name_gene_class_bottom{t} = num2str(c(i,:));
  end
end
t = t+1;
name_gene_class{t} = '0';
% name_gene_class_bottom{t} = '0';

% for i = 1:2:t
%   name_gene_class_top{i} = ' ';
% end
% for i = 2:2:t
%   name_gene_class_bottom{i} = ' ';
% end

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
setdiff(1:t,C)
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

sX = X(:,idx);
sA = A(:,idx);

figure(fig_nr+2), imagesc(sX), colormap(gray)
ticks = cumsum([1; a_counts(1:end-1)]);
xticks(ticks)
ticklabels = name_gene_class(C);
xticklabels(ticklabels)
xtickangle(90)

col = [0.25 0.4 0.7 0.9];
rep = 100;
imZclust = [col(1)*repmat(Zet(:,1), 1,rep) col(2)*repmat(Zet(:,2), 1,rep) ...
            col(3)*repmat(Zet(:,3), 1,rep) col(4)*repmat(Zet(:,4), 1,rep)];

        
rowA = unique(sA,'rows','stable'); % Finding replicates. 
rowA = rowA([1 2 3 4],:);
figure(27), imagesc(rowA), colormap(gray), title('Fisher unique')
      
imA = [col(4)*repmat(rowA(1,:),25,1); col(3)*repmat(rowA(2,:),25,1);...
       col(2)*repmat(rowA(3,:),25,1); col(1)*repmat(rowA(4,:),25,1)];
   
imA = [imA zeros(100,4*rep)];
imA = ones(size(imA))*0.1 + imA;

imZclust = ones(size(imZclust))*0.1 + imZclust;

figure(32), imagesc(imA)
figure(30), imagesc([imA; 0.6*sX imZclust]), colormap(jet(64))
xticks(ticks)
xticklabels(ticklabels)
xtickangle(90)
title('P3CL 1% filtering')
xlabel('Fisher hypothesis tests, p-value = 0.05, rejecting 13224 H0`s (93.3%)')









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



    