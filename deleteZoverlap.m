function Zredux = deleteZoverlap(Z, min_class, fig_nr)

% Finds and delets candidate columns that can be represented by a boolean
% sum of other candidate columns

if ~islogical(Z)
  disp('Logical, please')
  return
end

[n, nZ] = size(Z);
sumZ = sum(Z,1);
[~, idx] = sort(sumZ,'descend');
Z = Z(:,idx);

keep = false(1,nZ);
parfor k = 1 : nZ % Check the columns one by one
  Zk = Z;
  cand = true(1,nZ);
  % If the other column overlaps the 0's in the column in question, it is
  % not a representative
  for i = k+1 : nZ
    if sum(Zk(~Zk(:,k),i)) > min_class
      cand(i) = false;
    end
  end
  cand(1:k) = false;
  
  Zcand = Z(:,cand);
  repr = logical(sum(Zcand,2));
%   if fig_nr
%     subplot(1,2,1), imagesc([Z(:,k) repr (Z(:,k) & repr)]), colormap(gray), 
%     title('Original Representation Sum')
%   end
  if sum(Z(:,k) & repr) <  sum(Z(:,k))
    keep(k) = true;
  end
end
delete(gcp('nocreate'))
Zredux = Z(:,keep);

 

