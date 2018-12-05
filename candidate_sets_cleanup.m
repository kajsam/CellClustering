function [Zsets, sim_idx] = candidate_sets_cleanup(Zsets, X, Z, min_class, sim_idx)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 8th 2018

% Requires:     similar_cell.m, 

% The candidate sets need some clean-up. Sets of size 1 is deleted.
% Replicate sets are deleted. Sets that have none or more than one 1 for a
% cell are cleaned so that all cells have one cluster label. This takes
% time, but must be done. 

[n,d] = size(X);

% Order the sets
n_clust = zeros(1,size(Zsets,2));
for z = 1: size(Zsets,2)
  n_clust(z) = size(Zsets{z},2);
end
[n_clust,idx] = sort(n_clust,'descend');
Zsets = Zsets(idx);

% Sets of size 1
del = false(1, size(Zsets,2));
for z = 1: size(Zsets,2)
  if length(Zsets{z}) < 2
    del(z) = true;
  end
end
Zsets(del) = [];
n_clust(del) = [];

% Identical sets
del = false(1, size(Zsets,2));
for ii = min(n_clust):max(n_clust)
  idx = find(n_clust == ii);
  for z = 1: length(idx)
    for zz = z+1: length(idx)
      if isempty(setdiff(Zsets{idx(z)},Zsets{idx(zz)}))
        del(idx(z)) = true;
      end
    end
  end
end
Zsets(del) = [];

for z = 1: size(Zsets,2)
  Zsets{z} = Z(:,Zsets{z});
end

% if isempty(sim_idx)
%   tic
%   sim_idx = similar_cell(X);
%   toc
% end
% 
% for z = 1: size(Zsets,2)
%   Zsets{z} = label_cells(Z(:,Zsets{z}),X, sim_idx);
% end
% 
% del = false(1, size(Zsets,2));
% for z = 1: size(Zsets,2)
%   wons = sum(Zsets{z});
%   noise = find(wons<min_class); 
%   if ~isempty(noise)
%     Zsets{z}(:,noise) = [];
%     if size(Zsets{z},2)<2
%       del(z) = true;
%     end
%     Zsets{z} = label_cells(Zsets{z},X, sim_idx);
%   end
% end
% Zsets(del) = [];
