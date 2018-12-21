function Z = cond_prob_pi(X, pi, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), November 14th 2018



% This function calculates the conditional probability for each pair of
% rows, incorporating the gene probability.

% It has been a real pain in the ass: X(l,:).*pi(i,:) works, with or
% without the normalization. But X(l,:).*pi(i,:) is the likelihood of pi
% when X is observed to be 1, and this does not make sense, it should be
% the other way around: X(i,:).*pi(l,:), which doesn't work. 

% Input:    X - binary gene expression matrix
%           pi - gene probability

% Output:   Z - pairwise conditional probability

% Requires:      nothing

[n,d] = size(X);

if ~islogical(X)
  disp('Logical, please')
  return
end

% The conditional probability matrix
Z = zeros(n,n);

sumX = sum(X,2);
sumX(sumX==0) = 1;
parfor i = 1: n  
  XnX = zeros(1,n);
  % For each row, compare to all other rows 
  for l = 1 :n
    XnX(l) = sum(X(l,:).*pi(i,:));   
  end   
  
  Z(:,i) = XnX./sumX(i);        % Normalised by X(i) or by pi(i,:) or not
end

if fig_nr
  imZ = Z;
  imZ(logical(eye(n))) = median(Z);
  figure(fig_nr), imagesc(imZ), colormap(gray)
  title('Paired conditional success probability of cells')
  ylabel('Success in cell $$\ell$$ predicts success in cell $$i$$s','interpreter','latex')
  xlabel('Success in cell $$i$$ predicted by success in cell $$\ell$$s','interpreter','latex')
end 
% delete(gcp('nocreate'))

