function logZet = relabel_likelihood_cellcell(X, A, Zet, fig_nr)

% Input:    X - binary gene expression matrix
%           A - binary structure matrix approximation
%           cell_effect

if islogical(X) && islogical(A)
  if ~all(size(X)==size(A))
    disp('Wrong dimensions')
    return
  end
else
  disp('Logical, please')
  return
end

[n,d] = size(X);
K = size(Zet,2);
[~,idx] = sort(sum(Zet),'descend');
Zet = Zet(:,idx);

Au = false(K,d);
for k = 1: K
  id = find(Zet(:,k),1);
  Au(k,:) = A(id,:);
end

sum0 = sum(Au,2);
Au = Au(logical(sum0),:);

% Cell effect
XnA = X & A; 
cell_effect = sum(XnA,2); %(A)
w = cell_effect./sum(A,2); %d; %
m = (1-w); 
% figure(fig_nr+1), imagesc([w m]), title('Cell effect'), drawnow

% Gene effect
pi = zeros(K,d);
for k = 1: K
  gene_effect = sum(X(Zet(:,k),:)); %(Zet(:,k),:)
  gene_effect(gene_effect==0) = 1;
  gene_effect(gene_effect == sum(Zet(:,k))) = sum(Zet(:,k))-1;
  pi(k,:) = gene_effect/sum(Zet(:,k)); % n; %
end
% figure(fig_nr+2), imagesc(pi), title('Gene effect'), drawnow

% Maximum log likelihood
logL = zeros(n,K);
logZet = Zet;
for i = 1 :n
  for k = 1: K
    logL(i,k) = sum(log(w(i)*pi(k,X(i,:) & Au(k,:)))) ...
              + sum(log((1-w(i)*pi(k,~X(i,:) & Au(k,:)))))...% % sum(log((m(i))*(1-pi(k,~X(i,:) & Au(k,:)))))... % 
              + sum(log(w(i)*pi(k,X(i,:) & ~Au(k,:)))) ...
              + sum(log((1-w(i)*pi(k,~X(i,:) & ~Au(k,:))))); % sum(log((m(i))*(1-pi(k,~X(i,:) & ~Au(k,:))))); % 
  end
  [~,mix] = max(logL(i,:));
  logZet(i,:) = false;
  logZet(i,mix) = true;
  % Gene effect
  pi = zeros(K,d);
  for k = 1: K
    gene_effect = sum(X(Zet(:,k),:)); %(Zet(:,k),:)
    gene_effect(gene_effect==0) = 1;
    gene_effect(gene_effect == sum(Zet(:,k))) = sum(Zet(:,k))-1;
    pi(k,:) = gene_effect/sum(Zet(:,k)); % n; %
  end
end

if fig_nr
  figure(fig_nr), subplot(1,2,1), imagesc(logL), colormap(gray)
  title('loglik')
end 

logZet = false(n,K);
[~,mix] = max(logL,[],2);
for i = 1 :n
  logZet(i,mix(i)) = true;
end

sum(logZet)
logZet = logZet(:,logical(sum(logZet)));
sum(logZet)
if fig_nr  
  figure(fig_nr), subplot(1,2,2), imagesc(logZet) 
  colormap(gray),title(strcat(num2str(size(logZet,2)), ' candidate columns'))
end

