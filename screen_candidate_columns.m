function Z = screen_candidate_columns(X, Z)

% Kajsa Mollersen (kajsa.mollersen@uit.no), November 8th 2018

% This function aims at reducing the number of candidate cell cluster. A
% Fisher test is performed for each gene. Those candidate cell clusters
% that are never the best are deleted. 

% Input:        X - binary matrix
%               Z - the candidate columns

nZ = size(Z,2);   % The number of candidate columns

[n,d] = size(X);
Pset = zeros(nZ,d);
parfor z = 1: nZ
  p = zeros(1,d);
  for j = 1: d
    column = X(:,j);       
    contab = [sum(~column(Z(:,z))) sum(column(Z(:,z)));
                sum(~column(~Z(:,z))) sum(column(~Z(:,z)))];
    [~,p(j)] = fishertest(contab);     
  end
   
  Pset(z,:) = p;
end
delete(gcp('nocreate'))

[~,idx] = min(Pset);

Z = Z(:,unique(idx));
