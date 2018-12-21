function [Zet, cprototype] = seed_clusters(Z, cprototype, min_class, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), December 6th 2018

% Input:    Z - conditional probability for each cell
%           min_class - minimum number of cells in a cluster

% Output:   Zet - clusters
%           cprototype - a prototype cell for each cluster

% Requires: nothing     

% This is a clustering algorithm. Each cluster consists of a prototype and
% its members. 

[n,~] = size(Z); % Number of cells

if isempty(cprototype) % Start from scratch

  % Pick n/2 seeds (potential prototypes) at random, the rest will be lost
  % lambs waiting to be assigned
  seed = randsample(n,floor(n/2));
  lamb = setdiff(1:n,seed)';
  memb = assign(Z, seed, lamb); % Assign each lamb to a seed

  % Any seeds without lambs is redefined to be a lamb
  prototype = seed(unique(memb));
  lamb = setdiff(1:n, prototype);
  memb = assign(Z, prototype, lamb); % Assign each lamb to a prototype

  % A lamb can give equal conditional probability; need to clean that up
  while length(unique(memb)) < length(prototype) 
    prototype = prototype(unique(memb));
    lamb = setdiff(1:n, prototype);
    memb = assign(Z, prototype, lamb); % Assign each lamb to a prototype
  end

  %% Repeat the procedure where the initial seed-lamb assignment is reversed
  rseed = setdiff(1:n,seed)';
  rlamb = seed;
  memb = assign(Z, rseed, rlamb); % Assign each lamb to a seed

  rprototype = rseed(unique(memb));
  rlamb = setdiff(1:n, rprototype);
  memb = assign(Z, rprototype, rlamb); % Assign each lamb to a prototype

  while length(unique(memb)) < length(rprototype) 
    rprototype = rprototype(unique(memb));
    rlamb = setdiff(1:n, rprototype);
    memb = assign(Z, rprototype, rlamb); % Assign each lamb to a prototype
  end

  %% Combine the prototypes and reassign
  cprototype = sort([prototype; rprototype]);
end
clamb = setdiff(1:n, cprototype);
memb = assign(Z, cprototype, clamb); % Assign each lamb to a prototype

cprototype = cprototype(unique(memb));
clamb = setdiff(1:n, cprototype);
memb = assign(Z, cprototype, clamb); % Assign each lamb to a prototype

while length(unique(memb)) < length(cprototype) % assign again
  cprototype = cprototype(unique(memb));
  clamb = setdiff(1:n, cprototype);
  memb = assign(Z, cprototype, clamb); % Assign each lamb to a prototype
end

%% Deleting clusters of size smaller than min_class

[~,~,ic] = unique(memb); % Counting the members
a_counts = accumarray(ic,1);

for thresh = 1: min_class
  cprototype = cprototype(a_counts > thresh);
  clamb = setdiff(1:n, cprototype);
  memb = assign(Z, cprototype, clamb); % Assign each lamb to a prototype
  
  while length(unique(memb)) < length(cprototype) 
    cprototype = cprototype(unique(memb));
    clamb = setdiff(1:n, cprototype);
    memb = assign(Z, cprototype, clamb); % Assign each lamb to a prototype
  end
  [~,~,ic] = unique(memb);
  a_counts = accumarray(ic,1);
end

% Make the prototypes members
for i = 1: n
  memb(i) = find(Z(cprototype,i) == max(Z(cprototype,i)),1);
end

% Gather the cells as clusters
K = length(cprototype);
Zet = false(n,K);
for k = 1: K
  Zet(memb==k,k) = true;
end

if fig_nr  % Display the prototypes and cluster size
  table(a_counts,'RowNames',strsplit(num2str(cprototype')))
  imZ = Z(cprototype,cprototype);
  for i = 1: length(cprototype)
    imZ(i,i) = median(median(Z(cprototype,cprototype)));
  end
  figure(fig_nr), imagesc(imZ), colormap(jet)
  title('Conditional probability of prototypes')
end


function memb = assign(Z, seed, lamb)
% For each lamb, assign it to the seed that gives highest conditional 
% probability
memb = zeros(1, length(lamb)); % Each lamb is assigned a membership
for i = 1: length(lamb)
  memb(i) = find(Z(seed,lamb(i)) == max(Z(seed,lamb(i))),1);
end
