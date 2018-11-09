function [Zsets, sim_idx] = candidate_sets_cleanup(Zsets, X, Z, min_class, sim_idx)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 8th 2018

% The candidate sets need some clean-up. Sets of size 1 is deleted.
% Replicate sets are deleted. Sets that have none or more than one 1 for a
% cell are cleaned so that all cells have one cluster label. This takes
% time, but must be done. 

[n,d] = size(X);

size(Zsets,2)
del = false(1, size(Zsets,2));
for z = 1: size(Zsets,2)
  if length(Zsets{z}) < 2
    del(z) = true;
  end
end

Zsets(del) = [];
size(Zsets,2)

% Identical sets
del = false(1, size(Zsets,2));
for z = 1: size(Zsets,2)
  for zz = z+1: size(Zsets,2)
    if isempty(setdiff(Zsets{z},Zsets{zz}))
      del(z) = true;
    end
  end
end

Zsets(del) = [];
size(Zsets,2)

if sim_idx == 0
tic
sim_idx = zeros(1,n);
parfor ii = 1: n
  FPN = d*ones(1,n);
  for i = 1:n
    FPN(i) = sum(xor(X(ii,:),X(i,:)));                
  end
  FPN(ii) = d;
  [~,sim_idx(ii)] = min(FPN);
end
delete(gcp('nocreate'))
toc
end

for z = 1: size(Zsets,2)
  Zsets{z} = label_cells(Z(:,Zsets{z}),X, sim_idx);
  z
end

del = false(1, size(Zsets,2));
for z = 1: size(Zsets,2)
  wons = sum(Zsets{z});
  noise = find(wons<min_class);
  if ~isempty(noise)
    Zsets{z}(:,noise) = [];
    if size(Zsets{z},2)<2
      del(z) = true;
    end
    Zsets{z} = label_cells(Zsets{z},X, sim_idx);
  end
  % figure(fig_nr), imagesc(Zsets{z}), colormap(gray)
end
Zsets(del) = [];

size(Zsets,2)

