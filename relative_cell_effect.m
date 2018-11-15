function [A, idx] = relative_cell_effect(X, Zsets, A, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 12th 2018

% Input:        X - binary matrix
%               Zsets - the sets
%               A - fisher matrix

[n,d] = size(X);
K = size(Zsets,2)
[~,sidx] = sort(sum(Zsets),'descend');
Zsets = Zsets(:,sidx);
figure(fig_nr), imagesc(Zsets), colormap(gray), title('Classes')

cell_effect = sum(X,2);
% w = cell_effect/max(cell_effect);
% m = (1-cell_effect)/max(1-cell_effect);
w = cell_effect/d; 
m = (1-w); 
imZ = double(Zsets);
imZ2 = double(Zsets);

for k = 1: K
  Zet = Zsets(:,k);
  id = find(Zet,1);
  
  Aclass = A(id,:);

  Xclass = X(Zet,:);
  [nclass,d] = size(Xclass);

  % This is the gene effect
  gene_effect = sum(Xclass);
  gene_effect(gene_effect==0) = 1;
  gene_effect(gene_effect == nclass) = nclass-1;
  pi = gene_effect/nclass;
  figure(fig_nr+1), imagesc(repmat(gene_effect, 2,1)), colormap(gray)
  title(strcat('gene effect class ',num2str(k)))

  cell_effect = sum(Xclass,2);
  figure(fig_nr+2), imagesc(repmat(cell_effect, 1, 2)), colormap(gray)
  title(strcat('cell effect class',num2str(k)))

  % Log likelihood
  
  logL1 = zeros(1,nclass);
  logL0 = zeros(1,nclass);
  for i = 1: nclass
    logL1(i) = sum(log(w(i)*pi(Xclass(i,:) & Aclass))) ...
              + sum(log((m(i))*(1-pi(~Xclass(i,:) & Aclass))));
    logL0(i)         = sum(log(w(i)*pi(Xclass(i,:) & ~Aclass))) ...
              + sum(log((m(i))*(1-pi(~Xclass(i,:) & ~Aclass))));
  end
  logL = logL1+logL0;
  figure(fig_nr+3), subplot(1,3,1), hist(logL, 25), title(strcat('log lik ',num2str(k)))
  subplot(1,3,2), hist(logL1, 25), title(strcat('A = 1 '))
  subplot(1,3,3), hist(logL0, 25), title(strcat('A = 0 '))
  
  L = logL - min(logL);
  L = L./max(L);
  thresh = otsu_thresh(L,256)
  
  fig_nr = fig_nr+4;
  
  idx{k} = find(L<thresh);
  idZ = find(Zet);
    
  imZ(idZ(idx{k}),k) = 1.25;
  imZ2(idZ,k) = logL;
end
figure, imagesc(imZ)
figure, imagesc(imZ2)
