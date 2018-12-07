function [probZ, Zet, A, Zclust] = ...
    cell_seed_main005(X, probZ, Zet, A, Zclust, fig_nr)

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

if isempty(probZ)
     
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

if isempty(probZ)
  disp('Pairwise conditional success probability ... ')
  tic
  [~, probZ] = cond_prob_matrix(X, min_class, fig_nr);
  toc
  return
end

if isempty(Zet)
  disp('Cluster by seeding ... ')
  
  [Zet, prototype] = seed_clusters(probZ, min_class, fig_nr);
  
  return   
end

if isempty(A)
 
  tic
  [A,~,h] = fisher_set(X,Zet,fig_nr+1);
  toc
  figure(fig_nr), imagesc(Zet)
  ylabel(sum(h))
  title('Before maximum log likelihood')
  drawnow
  return
end

if isempty(Zclust)
  
  Ar = A; 
  maxR = 50
  Zloglik = Zet;
  rZet = Zet;
  for r = 1: maxR
    if fig_nr > 20
      fig_nr = 1;
    end
    fig_nr = fig_nr+2;
  
    rZet = relabel_likelihood_cellcell(X, Ar, rZet, 0); %
    [Ar,Psetr,h] = fisher_set(X, Zet,0);
    
    figure(fig_nr+2), imagesc(rZet)
    title('Maximum likelihood set'), xlabel(r), 
    ylabel(round(100*sum(h)/d))
    drawnow
    % figure(fig_nr+2), imagesc(Ar), colormap(gray)
    % title('Fisher matrix'), drawnow
    
    if isequal(Zloglik,rZet) || (r == maxR)
      [r sum(h)]   
      Zclust = rZet;
      return
    else
      Zloglik = rZet;
    end
  end
end

%[sX, Zet, sAr] = cellsort001_idx1(X, Zclust, Ar);

%[sort_idx, Zet, ticks, ticklabels] = sort_and_display(sX, Zet, sAr, Psetr, fig_nr);
