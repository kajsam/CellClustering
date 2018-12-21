function [A,p,h] = fisher_set(X, Zsets, alpha, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 8th 2018

% Input:        X - binary matrix
%               Zsets - the clusters
%               alpha - the significance of the hypothesis test

% Output:       A - structure matrix
%               p - the p-value for each gene
%               h - the reject or not for each gene

% Requires:     nothing


[n,d] = size(X);
% Decides if the gene is an all 0 or all 1 in case the hypotheses are
% rejected
thresh = mean(mean(X))*n; 

A = false(n,d);
p = zeros(1,d);
h = false(1,d);
parfor j = 1: d     % for each gene
  column = X(:,j);  % the observations for that gene
  Zset = Zsets;     % assigning for parallel computing
  K = size(Zset,2); % number of clusters
  
  % Order the clusters according to their mean value
  pi_class = zeros(1,K);
  for k = 1: K
    pi_class(k) = mean(column(Zset(:,k)));
  end
  [~,idx] = sort(pi_class);
  Zset = Zset(:,idx);
  
  % Do Fisher hypothesis test between the clusters
  pk = ones(1,K-1);
  hk = ones(1,K-1);
  for k = 1: K-1
    low = logical(sum(Zset(:,1:k),2)); % Concatenate those with lowest mean
    high = logical(sum(Zset(:,k+1:K),2)); % Those with highest mean
          
    contab = [sum(~column(low)) sum(column(low));
              sum(~column(high)) sum(column(high))];
    [hk(k),pk(k)] = fishertest(contab,'Alpha', alpha);     
  end
  % Find the clustering of clusters that gave the lowest p-value, and
  % assign to the structure matrix accordingly
  [minp,idx] = min(pk); 
  p(j) = pk(idx);
  h(j) = hk(idx);
  if minp < alpha 
    A(:,j) = logical(sum(Zset(:,idx+1:K),2));
  else % if the hypothesis was not rejected
    if sum(column) > thresh
      A(:,j) = true;
    end
  end  
end
if fig_nr
  figure(fig_nr), imagesc(A), colormap(gray), title(round(100*sum(h)/d))
  xlabel(mean(p))
  drawnow  
end
      