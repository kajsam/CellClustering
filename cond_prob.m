function Z = cond_prob(X)

% Kajsa Mollersen (kajsa.molllersen@uit.no), December 12th 2018

% This function calculates the conditional probability for each pair of
% rows. 

% Input:    X - binary gene expression matrix

% Output:   Z - the pairwise conditional probability

% Requires: nothing

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
    XnX(l) = sum(X(l,:)& X(i,:));      % Which entries are both 1
  end           
  Z(:,i) = XnX./sumX(i);        % Normalised by the total number of 1's in the row 
end
% delete(gcp('nocreate'))
