function Zsets = candidate_sets(Z, delete_noise, min_class, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no), November 8th 2018

% Creates a set of columns from Z for each single column of Z. The set will
% ideally have a single 1 for each cell, but this is normally not achieved.
% A clean-up follows.

% Input:    Z : candidate columns

% Output:   Zsets : indexes for each set

if ~islogical(Z)
  disp('Logical, please')
  return
end

[n, nZ] = size(Z);
Zsets = cell(1,nZ);

for k =  1: nZ   
  w = Z(:,k);                           % for each column
  set = [];                             % create a set
  for r = 1 : 50                       % maximum size of the set
    w_or_set = logical(sum([w Z(:,set)],2));     % row sum of the set       
    
    pen = ones(1,nZ)*n;                 % penalties
    rew = zeros(1,nZ);                  % rewards
    crit = ones(1,nZ)*(-n);             % selection criterion
    
    for i = setdiff(1:nZ,[k set])                   % for columns not in the set
      pen(i) = sum((w_or_set + Z(:,i))>1);          % penalize for row sum
      rew(i) = sum(~logical(w_or_set) & Z(:,i));    % reward for previous 0 entry
      crit(i) = rew(i) - pen(i);                    % balance penalty and reward   
    end
    
    if ~isempty(find(pen== 0,1))        % if a column produces zero penalty
      crit(pen > 0) = -n;                % all others are of no interest
    end
    [maxcrit, idx] = max(crit);      
    
    if maxcrit < 0                      % now it is time to stop
      Zet = [k set];                % not including the lates column

      if delete_noise                   % delete noise columns
        del = false(1,length(Zsets{k}));
        for i = 1: length(Zsets{k})
          if sum(Z(:,Zsets{k}(i))) - sum(Z(:,Zsets{k}(i)) & logical(sum(Z(:,setdiff(Zsets{k},Zsets{k}(i))),2))) < min_class
            del(i) = true;
          end
        end
        Zet(del) = [];
      end
      
      Zsets{k} = Zet;
      
      if fig_nr
        figure(fig_nr), imagesc(Z(:,Zsets{k})), colormap(gray), xlabel(k)
        title(k), drawnow
      end
      
      break                             % stop if maxcrit < 0
      
    else
      set = [set idx];                  % add to set
    end
    if r == 50
      disp('r = 50!')
      return
    end
  end
end