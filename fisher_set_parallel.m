function Pset = fisher_set_parallel(X, Zsets)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 8th 2018

% Input:        X - binary matrix
%               Zsets - the sets

nZets = size(Zsets,2);   % The number of sets we are dealing with

[n,d] = size(X);
alpha = 0.05;
Pset = zeros(nZets,d);
parfor z = 1: nZets
  Zset = Zsets{z};
  K = size(Zset,2);
  pi_class = zeros(1,K);
  p = zeros(1,d);
  for j = 1: d
    column = X(:,j);
    for k = 1: K
      pi_class(k) = mean(column(Zset(:,k)));
    end
    [pi_class,idx] = sort(pi_class);
    Zset = Zset(:,idx);
    pk = ones(1,K-1);
    for k = 1: K-1
      low = logical(sum(Zset(:,1:k),2));
      high = logical(sum(Zset(:,k+1:K),2));
          
      contab = [sum(~column(low)) sum(column(low));
                sum(~column(high)) sum(column(high))];
      [~,pk(k)] = fishertest(contab,'Alpha', alpha);     
    end
    [~,idx] = min(pk);
    p(j) = pk(idx);
  end
  Pset(z,:) = p;
end
delete(gcp('nocreate'))
      