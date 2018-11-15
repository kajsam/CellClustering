function [binZ, thresh] = cond_prob_matrix(X, min_class, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), November 14th 2018

% Requires:      otsu_thresh.m

% This function calculates the conditional probability for each pair of
% rows. Then, the Otsu-threshold is calculated for each row, and a
% binarization is done. The result is a set of columns that serves as
% candidate columns for the cell clustering.


% Input:    X - binary gene expression matrix
%           min_class - minimum size of cell cluster

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
  Xi = X(i,:);
  XnX = zeros(1,n);
  % For each row, compare to all other rows 
  for l = 1 :n
    XnX(l) = sum(X(l,:)& Xi);      % Which entries are both 1
  end           
  Z(:,i) = XnX./sumX(i);        % Normalised by the total number of 1's in the row 
end
delete(gcp('nocreate'))

if fig_nr
  imZ = Z;
  imZ(logical(eye(n))) = median(Z);
  figure(fig_nr), subplot(1,2,1), imagesc(imZ), colormap(gray)
  title('Paired conditional success probability of cells')
  ylabel('Success in cell $$\ell$$ predicts success in cell $$i$$s','interpreter','latex')
  xlabel('Success in cell $$i$$ predicted by success in cell $$\ell$$s','interpreter','latex')
end 

% Calculate the threshold for binarization. The diagonal is excluded.  
num_bins = 256; 
thresh = zeros(1,n);
binZ = true(n,n);

for i = 1: n
  z = Z(i,[1:i-1 i+1:n]) - min(Z(i,[1:i-1 i+1:n]));
  if max(z) ~= 0
    z = z./max(z);
  end
  thresh(i) = otsu_thresh(z,num_bins);
  z(z<thresh(i)) = 0;
  binZ(i,[1:i-1 i+1:n]) = logical(z);
end

binZ = unique(binZ', 'rows', 'stable'); % Finding replicates. 
binZ = binZ';

% Delete the clusters that are smaller than the smallest class size
sumZ = sum(binZ,1);
idx0 = sumZ < min_class;    
binZ(:,idx0) = [];               

if fig_nr  
  figure(fig_nr), subplot(1,2,2), imagesc(binZ) 
  colormap(gray), title('Binarised by row'), 
  xlabel('Each column represents a proposed cluster')
  drawnow
end
