function [Z, redZ, Zsets, sim_idx, Pset, Zclust, Ar, Psetr, sort_idx] = ...
    cell_clustering_main001(X, Z, redZ, Zsets, sim_idx, Pset, Zclust, Ar, Psetr, sort_idx, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), November 8th 2018

% This is the main file for a method that clusters cell in a binary
% single-cell RNAseq matrix

% Requires:     cond_prob_matrix.m, otsu_thresh.m, candidate_sets.m, 
%               candidate_sets_cleanup.m, similar_cell.m, label_cells.m, 
%               fisher_set_parallel.m, fisher_set.m, relabel_likelihood.m,
%               sort_and_display.m, relative_cell_effect.m

% Input:        X - binary matrix of gene expression
%                  
% Optional:     Z0 - candidate columns
%               redZ0 - reduced set of candidate columns
%               Zsets0 - column sets
%               sim_idx - most similar cell for each of the cells in X
%               Pset0 - p-values for the column sets

[n,d] = size(X);

min_class = ceil(0.01*n); % Minimum class size. 

if isempty(Z)
     
  cell_effect = sum(X,2);
  cell_effect = cell_effect/d;
  
  gene_effect = sum(X,1);
  gene_effect = gene_effect/n;
  
  imcell = cell_effect-min(cell_effect);
  imcell = imcell/max(imcell);
  imX = [X repmat(imcell,1,50)];
  imgene = [gene_effect zeros(1,50)];
  imX = [imX; repmat(imgene,50,1)];
  figure(fig_nr), imagesc(imX), colormap(gray)
  title('Observed values')
  % set(gca,'xaxisLocation','top')
  ylabel('Cells')
  xlabel('Genes')
  fig_nr = fig_nr+1;
end

%% First time is without cell and gene effect

if isempty(Z)
  disp('Pairwise conditional success probability ... ')
  tic
  Z = cond_prob_matrix(X, min_class, fig_nr);
  toc
  return
end

if isempty(redZ)
  disp('Reducing the number of cluster proposals ... ')
  
  % This can be done to save time

  % redZ = deleteZoverlap(Z, min_class, 0); 
  % figure(fig_nr+1), imagesc(redZ), colormap(gray)
  % title('Reduced set'), drawnow

  % redZ = screen_candidate_columns(X, redZ);
  % figure(fig_nr+1), imagesc(redZ), colormap(gray)
  % title('Screened set'), drawnow
  
  return   
end

sumZ = sum(redZ,1);
[~, idx] = sort(sumZ,'descend');
redZ = redZ(:,idx);

if isempty(Zsets)
  disp('Creating sets of clusters')
  tic
  Zsets = candidate_sets(redZ, 1, min_class, 0);
  toc
%   % This is a quality measure for the column sets. Can save a lot of time
%   qual = zeros(1, length(Zsets));
%   for s = 1: length(Zsets)
%     q = sum(redZ(:,Zsets{s}),2);
%     qual(s) = sum(~q)+ sum(q>1);
%   end
%   [~, qidx] = sort(qual,'ascend');
%   Zsets_sort = Zsets(qidx);
%   Zsets_half = Zsets_sort(1:floor(length(Zsets_sort)));
  disp('Single labeling ... ')
  tic
  [Zsets, sim_idx] = candidate_sets_cleanup(Zsets, X, redZ, min_class, sim_idx);
  toc
  return
end

if isempty(Pset)
  disp('Fisher matrices (takes time) ...')
  tic
  Pset = fisher_set_parallel(X, Zsets);
  toc
  return
end

if isempty(Zclust)
  p = sum(Pset,2);
  [~, idx] = min(p)
  
  for idx = 29
      Zet = Zsets{idx};
  figure(fig_nr), imagesc(Zet), title(idx)
  pause
  end
 
  [A,~,h] = fisher_set(X,Zet,fig_nr+1);
  Ar = A; 
  maxR = 50
  Zloglik = Zet;
  for r = 1: maxR
    if fig_nr > 20
      fig_nr = 1;
    end
    fig_nr = fig_nr+1;
  
    Zet = relabel_likelihood(X, Ar, Zet, 0); %
    
    figure(fig_nr), imagesc(Zet), colormap(gray)
    title('Maximum likelihood set'), xlabel(r), drawnow
    
    [Ar,Psetr,h] = fisher_set(X,Zet,0);
    % figure(fig_nr+2), imagesc(Ar), colormap(gray)
    % title('Fisher matrix'), drawnow
    
    if isequal(Zloglik,Zet) || (r == maxR)
        r
        sum(h)
        [sort_idx, Zet] = sort_and_display(X, Zet, Ar, Psetr, fig_nr);
        Zclust = Zet;
      return
    else
      Zloglik = Zet;
    end
  end
end

% if isempty(sort_idx)

[sX, Zet, sAr] = cellsort001_idx29(X, Zclust, Ar);

[sort_idx, Zet, ticks, ticklabels] = sort_and_display(sX, Zet, sAr, Psetr, fig_nr);

%   return
%end


X = X(:,sort_idx);
sAr = Ar(:,sort_idx);








