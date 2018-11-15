function [A,Pset,h] = fisher_set(X, Zsets, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 8th 2018

% Input:        X - binary matrix
%               Zsets - the sets

if iscell(Zsets)
  nZets = size(Zsets,2);   % The number of sets we are dealing with
else
  nZets = 1;
  Zsets0 = Zsets;
  Zsets = cell(1,1);
  Zsets{1} = Zsets0;
end
[n,d] = size(X);
thresh = mean(mean(X))*n;
alpha = 0.05;
Pset = zeros(nZets,d);
for z = 1: nZets
  Zset = Zsets{z};
  K = size(Zset,2);
  pi_class = zeros(1,K);
  A = false(n,d);
  p = zeros(1,d);
  h = false(1,d);
  for j = 1: d
    column = X(:,j);
    for k = 1: K
      pi_class(k) = mean(column(Zset(:,k)));
    end
    [pi_class,idx] = sort(pi_class);
    Zset = Zset(:,idx);
    pk = ones(1,K-1);
    hk = ones(1,K-1);
    for k = 1: K-1
      low = logical(sum(Zset(:,1:k),2));
      high = logical(sum(Zset(:,k+1:K),2));
          
      contab = [sum(~column(low)) sum(column(low));
                sum(~column(high)) sum(column(high))];
      [hk(k),pk(k)] = fishertest(contab,'Alpha', alpha);     
    end
    [minp,idx] = min(pk);
    p(j) = pk(idx);
    h(j) = hk(idx);
    if minp < alpha
      A(:,j) = logical(sum(Zset(:,idx+1:K),2));
    else
      if sum(column) > thresh
        A(:,j) = true;
      end
    end  
  end
  if fig_nr
    figure(fig_nr), imagesc(Zset), colormap(gray), title('Candidate set')
    xlabel(z)
    figure(fig_nr+1), imagesc(A), colormap(gray), title(sum(h))
    xlabel(mean(p))
    drawnow  
  end
  Pset(z,:) = p;
end
      