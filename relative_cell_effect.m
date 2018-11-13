function [A] = relative_cell_effect(X, Zsets, A, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) November 12th 2018

% Input:        X - binary matrix
%               Zsets - the sets
%               A - fisher matrix

[n,d] = size(X);
K = size(Zsets,2)
[~,idx] = sort(sum(Zsets),'descend');
Zsets = Zsets(:,idx);
figure(fig_nr), imagesc(Zsets), colormap(gray), title('Classes')

% A(:,sum(A)==n) = false;
cell_effect = sum(X,2);
w = cell_effect/max(cell_effect);
m = (1-cell_effect)/max(1-cell_effect);
for k = 1: K
  Zet = Zsets(:,k);
  id = find(Zet,1);
  
  Aclass = A(id,:);

  Xclass = X(Zet,:);
  [nclass,d] = size(Xclass);

  % This is the gene effect
  gene_effect = sum(Xclass);
  pi = gene_effect/nclass;
  figure(fig_nr+1), imagesc(repmat(gene_effect, 2,1)), colormap(gray)
  title(strcat('gene effect class ',num2str(k)))

  cell_effect = sum(Xclass,2);
  %w = cell_effect/max(cell_effect);
  figure(fig_nr+2), imagesc(repmat(cell_effect, 1, 2)), colormap(gray)
  title(strcat('cell effect class',num2str(k)))

  % Archetype cell
  
  logL1 = zeros(1,nclass);
  logL0 = zeros(1,nclass);
  for i = 1: nclass
    logL1(i) = sum(log(w(i)*pi(Xclass(i,:) & Aclass))) ...
              + sum(log((m(i))*(1-pi(~Xclass(i,:) & Aclass))));
    logL0(i)         = sum(log(w(i)*pi(Xclass(i,:) & ~Aclass))) ...
              + sum(log((m(i))*(1-pi(~Xclass(i,:) & ~Aclass))));
  end
  logL = logL1+logL0;
  % [maxL,arch] = max(logL);
  figure(fig_nr+4), subplot(1,3,1), hist(logL, 25), title(strcat('log lik ',num2str(k)))
  subplot(1,3,2), hist(logL1, 25), title(strcat('A = 1 '))
  subplot(1,3,3), hist(logL0, 25), title(strcat('A = 0 '))

  
%   logL1 = zeros(1,n);
%   logL0 = zeros(1,n);
%   for i = 1: n
%     logL1(i) = sum(log(w(i)*pi(X(i,:) & Aclass))) ...
%               + sum(log((1-w(i))*(1-pi(~X(i,:) & Aclass))));
%     logL0(i)         = sum(log(w(i)*pi(X(i,:) & ~Aclass))) ...
%               + sum(log((1-w(i))*(1-pi(~X(i,:) & ~Aclass))));
%   end
%   logL = logL1+logL0;
%   [maxL,arch] = max(logL);
%   figure(fig_nr+5), subplot(1,3,1), hist(logL, 25), title(strcat('log lik ',num2str(k)))
%   subplot(1,3,2), hist(logL1, 25), title(strcat('A = 1 '))
%   subplot(1,3,3), hist(logL0, 25), title(strcat('A = 0 '))




%   warch = cell_effect(arch)/sum(Xclass(arch,:) & Aclass);    
%   class_val = zeros(1,nclass);
%   for i = 1: nclass
%     class_val(i) = sum(Xclass(i,:) & Aclass)*warch/cell_effect(i);    
%   end

  % figure(fig_nr+3), imagesc(repmat(class_val,2,1))
  % title('Class 1 validity')
%   figure(fig_nr+3), subplot(1,2,1), hist(class_val, 25), title(strcat('Class validity',num2str(k)))
%   warch = cell_effect(arch)/sum(~Xclass(arch,:) & ~Aclass);    
%   class_val = zeros(1,nclass);
%   for i = 1: nclass
%     class_val(i) = sum(~Xclass(i,:) & ~Aclass)*warch/cell_effect(i);    
%   end
%   figure(fig_nr+3), subplot(1,2,2), hist(class_val, 25), title(strcat('Class validity',num2str(k)))
  fig_nr = fig_nr+4;
  
  if k ==1
    idx = find(logL<-1500);
    idZ = find(Zet);
    imZ = double(Zet);
    imZ(idZ(idx)) = 1.25;
    
  elseif k ==2
    idx = find(logL<-1500);
    idZ = find(Zet);
    % imZ = double(Zet);
    imZ(idZ(idx)) = 0.5;
    
    figure, imagesc(imZ)
  end
end