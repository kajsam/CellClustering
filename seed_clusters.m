function [Zet, cprototype] = seed_clusters(Z, min_class, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), December 6th 2018

% Requires:      



% Input:    Z - conditional probability for each cell
[n,~] = size(Z);

% Pick n/2 seed cells at random, the rest will be lost lambs
seed = randsample(n,floor(n/2));
lamb = setdiff(1:n,seed);

% For each lost lamb, assign it to the seed that gives highest conditional 
% probability
memb = zeros(1, length(lamb)); % Each lost lamb becomes a member 
for i = 1: length(lamb)
  memb(i) = find(Z(seed,lamb(i)) == max(Z(seed,lamb(i))),1);
end

% Any seeds without lambs is redefined to be a lost lamb
prototype = seed(unique(memb));
prototype = sort(prototype)';
lamb = setdiff(1:n, prototype);
% For each lost lamb, assign it to the seed that gives highest conditional 
% probability
memb = zeros(1, length(lamb)); % Each lost lamb becomes a member 
for i = 1: length(lamb)
  memb(i) = find(Z(prototype,lamb(i)) == max(Z(prototype,lamb(i))),1);
end

if length(unique(memb)) < length(prototype) % assign again
  prototype = prototype(unique(memb));
  lamb = setdiff(1:n, prototype);
  memb = zeros(1, length(lamb)); % Each lost lamb becomes a member 
  for i = 1: length(lamb)
    memb(i) = find(Z(prototype,lamb(i)) == max(Z(prototype,lamb(i))),1);
  end
end

[C,~,ic] = unique(memb);
a_counts = accumarray(ic,1);
T = table(a_counts,'RowNames',strsplit(num2str(prototype)))

imZ = Z(prototype,prototype);
for i = 1: length(prototype)
  imZ(i,i) = median(median(Z(prototype,prototype)));
end
figure(fig_nr), imagesc(imZ), colormap(jet)

maxProb = zeros(1,length(prototype));
for i = 1: length(prototype)
  maxProb(i) = max(Z(prototype(i), lamb(memb ==i)));
end

% Reverse:
rseed = setdiff(1:n,seed);
rlamb = seed;
% For each lost lamb, assign it to the seed that gives highest conditional 
% probability

memb = zeros(1, length(rlamb)); % Each lost lamb becomes a member 
for i = 1: length(rlamb)
  memb(i) = find(Z(rseed,rlamb(i)) == max(Z(rseed,rlamb(i))),1);
end

% Any seeds without lambs is redefined to be a lost lamb
rprototype = rseed(unique(memb));
rprototype = sort(rprototype);
rlamb = setdiff(1:n, rprototype);
% For each lost lamb, assign it to the seed that gives highest conditional 
% probability
memb = zeros(1, length(rlamb)); % Each lost lamb becomes a member 
for i = 1: length(rlamb)
  memb(i) = find(Z(rprototype,rlamb(i)) == max(Z(rprototype,rlamb(i))),1);
end
if length(unique(memb)) < length(rprototype) % assign again
  rprototype = rprototype(unique(memb));
  rlamb = setdiff(1:n, rprototype);
  memb = zeros(1, length(rlamb)); % Each lost lamb becomes a member 
  for i = 1: length(rlamb)
    memb(i) = find(Z(rprototype,rlamb(i)) == max(Z(rprototype,rlamb(i))),1);
  end
end

[C,~,ic] = unique(memb);
a_counts = accumarray(ic,1);
T = table(a_counts,'RowNames',strsplit(num2str(rprototype)))

imZ = Z(rprototype,rprototype);
for i = 1: length(rprototype)
  imZ(i,i) = median(median(Z(rprototype,rprototype)));
end
figure(fig_nr), imagesc(imZ), colormap(jet)

maxProb = zeros(1,length(rprototype));
for i = 1: length(rprototype)
  maxProb(i) = max(Z(rprototype(i), rlamb(memb ==i)));
end


% Combine:
cprototype = sort([prototype rprototype]);
clamb = setdiff(1:n, cprototype);
% For each lost lamb, assign it to the seed that gives highest conditional 
% probability
memb = zeros(1, length(clamb)); % Each lost lamb becomes a member 
for i = 1: length(clamb)
  memb(i) = find(Z(cprototype,clamb(i)) == max(Z(cprototype,clamb(i))),1);
end
if length(unique(memb)) < length(cprototype) % assign again
  cprototype = cprototype(unique(memb));
  clamb = setdiff(1:n, cprototype);
  memb = zeros(1, length(clamb)); % Each lost lamb becomes a member 
  for i = 1: length(clamb)
    memb(i) = find(Z(cprototype,clamb(i)) == max(Z(cprototype,clamb(i))),1);
  end
end

[C,~,ic] = unique(memb);
a_counts = accumarray(ic,1);

[sa_counts, idx] = sort(a_counts,'descend');
sprototype = cprototype(idx);
T = table(sa_counts,'RowNames',strsplit(num2str(sprototype)))

imZ = Z(cprototype,cprototype);
for i = 1: length(cprototype)
  imZ(i,i) = median(median(Z(cprototype,cprototype)));
end
figure(fig_nr), imagesc(imZ), colormap(jet)

% Deleting clusters of size smaller than min_class
for thresh = 1: min_class
  cprototype = cprototype(a_counts > thresh);
  clamb = setdiff(1:n, cprototype);
  % For each lost lamb, assign it to the seed that gives highest conditional 
  % probability
  memb = zeros(1, length(clamb)); % Each lost lamb becomes a member 
  for i = 1: length(clamb)
    memb(i) = find(Z(cprototype,clamb(i)) == max(Z(cprototype,clamb(i))),1);
  end
  if length(unique(memb)) < length(cprototype) % assign again
    cprototype = cprototype(unique(memb));
    clamb = setdiff(1:n, cprototype);
    memb = zeros(1, length(clamb)); % Each lost lamb becomes a member 
    for i = 1: length(clamb)
      memb(i) = find(Z(cprototype,clamb(i)) == max(Z(cprototype,clamb(i))),1);
    end
  end

  [C,~,ic] = unique(memb);
  a_counts = accumarray(ic,1);
  T = table(a_counts,'RowNames',strsplit(num2str(cprototype)))

  imZ = Z(cprototype,cprototype);
  for i = 1: length(cprototype)
    imZ(i,i) = median(median(Z(cprototype,cprototype)));
  end
  figure(fig_nr), imagesc(imZ), colormap(jet)
end

for i = 1: n
  memb(i) = find(Z(cprototype,i) == max(Z(cprototype,i)),1);
end

K = length(cprototype);
Zet = false(n,K);
for k = 1: K
  Zet(memb==k,k) = true;
end

figure(fig_nr+1), imagesc(Zet)





