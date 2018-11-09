function [binZ, thresh] = association_matrix(X, min_class, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), November 8th 2018

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

% The association matrix
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
  Z(:,i) = XnX./sumX(i);        % Sum up row-wise, normalised by the total number of 1's in the row 
end
delete(gcp('nocreate'))
if fig_nr
  figure(fig_nr), subplot(1,2,1), imagesc(Z), colormap(gray)
  title('Association matrix')
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

sumZ = sum(binZ,1);

idx0 = sumZ < min_class;    % Smaller than the smallest class size
binZ(:,idx0) = [];                % Noise, delete it

if fig_nr  
  figure(fig_nr), subplot(1,2,2), imagesc(binZ) 
  colormap(gray), title(median(thresh)), drawnow
end
