function W = label_cells(W, X, sim_idx)

[n,d] = size(X);

n_sim = length(sim_idx);

for ii = 1: n_sim
  for t = 1: 5
    sumW = sum(W,2);
    k0 = find(sumW == 0);
    if isempty(k0)
      ii = n_sim;
      break
    end
    nk0 = length(k0);
    for k = 1: nk0  
      W(k0(k),W(sim_idx{ii}(k0(k)),:)~=0) = true;
    end
  end
end
  
sumW = sum(W,2);
k0 = find(sumW == 0);
if ~isempty(k0)
  nk0 = length(k0)
  idx0 = zeros(1, nk0);
  for k = 1: nk0
    FPN = d*ones(1,n);
    for i = setdiff(1:n, k0)
      FPN(i) = sum(xor(X(k0(k),:),X(i,:)));    
    end
    [~,idx0(k)] = min(FPN);
  end

  for k = 1: nk0  
    W(k0(k),W(idx0(k),:)~=0) = true;
  end
end

for ii = 1: n_sim
  for t = 1: 5
    sumW = sum(W,2);
    k1 = find(sumW > 1);
    if isempty(k1)
      ii = n_sim;
      break
    end
    nk1 = length(k1);
    for k = 1: nk1  
      W(k1(k),:) = W(sim_idx{ii}(k1(k)),:);
    end
  end
end

sumW = sum(W,2);
k1 = find(sumW > 1);

if ~isempty(k1)
  nk1 = length(k1)
  idx1 = zeros(1, nk1);

  for k = 1: nk1
    FPN = d*ones(1,n);
    for i = setdiff(1:n, k1)
      FPN(i) = sum(xor(X(k1(k),:),X(i,:)));    
    end
    [~,idx1(k)] = min(FPN);
  end

  for k = 1: nk1  
    W(k1(k),:) = W(idx1(k),:);
  end
end