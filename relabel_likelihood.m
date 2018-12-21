function [logZet, pi] = relabel_likelihood(X, A, Zet, min_class, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), December 20th 2018

% Input:    X - binary gene expression matrix
%           A - binary structure matrix 
%           Zet - clusters
%           min_class - minimum number of cells in a cluster

% Output:   logZet - updated clusters
%           pi - gene probability matrx

% Requires: nothing

if islogical(X) && islogical(A)
  if ~all(size(X)==size(A))
    disp('Wrong dimensions')
    return
  end
else
  disp('Logical, please')
  return
end

[n,d] = size(X);    % Number of cells and genes
K = size(Zet,2);    % Number of clusters

% Cell probability: The average number of genes expressed in each cell out
% of those genes that should be expressed according to the structure matrix
XnA = X & A; 
p_cell = sum(XnA,2)./sum(A,2); 

% Gene probability: The average number of cells expressed in a cluster in 
% each gene
pi_gene = zeros(K,d);
for k = 1: K
  gene_effect = sum(X(Zet(:,k),:)); % Sum of expressed cells for cluster k
  gene_effect(gene_effect==0) = 1; % If no cells are expressed
  % Or of all cellls are expressed
  gene_effect(gene_effect == sum(Zet(:,k))) = sum(Zet(:,k))-1;  
  pi_gene(k,:) = gene_effect/sum(Zet(:,k)); % By the cluster size
end

% Identify which genes should be expressed for each cluster
Au = false(K,d);
for k = 1: K
  id = find(Zet(:,k),1);
  Au(k,:) = A(id,:);
end

% Log likelihood for each cluster, assuming pi = pi_cell*pi_gene
logL = zeros(n,K);
for i = 1 :n
  for k = 1: K
    logL(i,k) = sum(log(p_cell(i)*pi_gene(k,X(i,:) & Au(k,:)))) ...
              + sum(log((1-p_cell(i)*pi_gene(k,~X(i,:) & Au(k,:)))))...
              + sum(log(p_cell(i)*pi_gene(k,X(i,:) & ~Au(k,:)))) ...
              + sum(log((1-p_cell(i)*pi_gene(k,~X(i,:) & ~Au(k,:)))));
  end

end

% Relabel cells according to the maximum log likelihood
logZet = false(n,K);
[~,mix] = max(logL,[],2);
for i = 1 :n
  logZet(i,mix(i)) = true;
end

% Discard clusters of size smaller than min_class
logL(:,sum(logZet)< min_class) = [];
% Relabel cells according to the maximum log likelihood
logZet = false(n,K);
[~,mix] = max(logL,[],2);
for i = 1 :n
  logZet(i,mix(i)) = true;
end
logZet = logZet(:,logical(sum(logZet)));

if fig_nr  
  figure(fig_nr), subplot(1,3,2), imagesc(logZet) 
  title('Loglik')
  ylabel(sum(logZet))
  set(get(gca,'YLabel'),'Rotation',0)
end
K = size(logZet,2);

% Recalculate the gene probability
pi = zeros(n,d);
for k = 1: K
  gene_effect = sum(X(logZet(:,k),:));
  gene_effect(gene_effect==0) = 1; 
  gene_effect(gene_effect == sum(Zet(:,k))) = sum(Zet(:,k))-1; 
  pi(logZet(:,k),:) = repmat(gene_effect/sum(logZet(:,k)),sum(logZet(:,k)),1);
end
