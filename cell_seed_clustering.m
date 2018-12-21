function [probZ, tZet, rA, prototypes] = ...
    cell_seed_clustering(X, probZ, prototypes, fig_nr)

% Kajsa Mollersen (kajsa.molllersen@uit.no), December 20th 2018

% This is the main file for a method that clusters cell in a binary
% single-cell RNAseq matrix

% Requires:     cond_prob.m, seed_clusters.m, fisher_set.m, relabel_likelihood.m,
%               cond_prob_pi.m, 

% Input:        X - binary matrix of gene expression
%                  
% Optional:     probZ - pairwise conditional probability of success
%               prototypes - a prototype cell for each cluster        
%
% Output:       probZ - pairwise conditional probability of success
%               tZet - the clusters
%               rA - structure matrix
%               prototypes - a prototype cell for each cluster       

tZet = []; % Empty output if not updated
rA = [];

[n,d] = size(X);    % n cells, d genes

min_class = ceil(0.01*n) % Minimum cluster size. 
alpha = 0.05 % For the Fisher hypothesis testing
maxR = 20 % Maximum number of iterations

% Displaying the observed values, and the mean cell and gene values   
mean_cell = sum(X,2)/d; 
mean_gene = sum(X,1)/n;
  
imcell = mean_cell-min(mean_cell);
imcell = imcell/max(imcell);
imX = [X repmat(imcell,1,round(d/50))];
imgene = [mean_gene zeros(1,round(d/50))];
imX = [imX; repmat(imgene,round(n/50),1)];
figure(fig_nr), imagesc(imX), colormap(gray)
title('Observed values')
ylabel('Cells')
xlabel('Genes')

if isempty(prototypes)
  
    %% First time is without cell and gene effect

    disp('Pairwise conditional success probability ... ')
    probZ = cond_prob(X);
 
    % Display
    imZ = probZ;
    imZ(logical(eye(n))) = median(probZ);
    figure(fig_nr+1), imagesc(imZ), colormap(gray)
    title('Paired conditional success probability of cells')
    ylabel('Success in cell $$\ell$$ predicts success in cell $$i$$s','interpreter','latex')
    xlabel('Success in cell $$i$$ predicted by success in cell $$\ell$$s','interpreter','latex')
    drawnow
    
    disp('Cluster by conditional probability ... ')
    [Zet, prototype] = seed_clusters(probZ, prototypes, min_class, 0);
    
    figure(fig_nr+2), imagesc(Zet), title('Very first clustering')
  
    disp('Fisher hypothesis tests ...')
    [A,~,~] = fisher_set(X,Zet,alpha, fig_nr+3);
        
    Ar = A; 
    Zloglik = Zet;
    rZet = Zet;
    protr = prototype;
    
    fig_nr = fig_nr+4;
    for r = 1: maxR
      if fig_nr > 20
        fig_nr = 4;
      end
      
      % Display the clustering
      figure(fig_nr), subplot(1,3,1), imagesc(rZet), 
      xlabel(r-1)
      
      % Calculate probability matrix
      [rZet, pi] = relabel_likelihood(X, Ar, rZet, min_class, fig_nr); %
      
      figure(fig_nr+1), imagesc(pi), colormap(gray)
      title('Gene probabilities')
    
      % Calculate conditional probability matrix
      piZ = cond_prob_pi_matrix(X, pi, 0);
    
      % Cluster by conditional probability
      [sZet, protr] = seed_clusters(piZ, protr, min_class, 0);
      
      figure(fig_nr), subplot(1,3,3), imagesc(sZet), 
      xlabel(r)
      drawnow
    
      % New Fisher matrix
      [Ar,~,~] = fisher_set(X, rZet,alpha, fig_nr+2);
      
      % Break if identical, or if you're bored
      if isequal(Zloglik,rZet) || (r == maxR)
        prototypes = protr;
        figure, imagesc(rZet), xlabel(r)
        break
      else
        Zloglik = rZet;
      end
      fig_nr = fig_nr+3;
      pause
    end
  return
end

% Repeat everything from above
maxT = 10;
protr = prototypes;
for t = 1: maxT
  [rZet, protr] = seed_clusters(probZ, protr, min_class, fig_nr);
  [rA,~,~] = fisher_set(X,rZet,alpha, 0);

  Zloglik = rZet;
  
  for r = 1 : maxR
    [rZet, pi] = relabel_likelihood(X, rA, rZet, min_class, 0); %
    piZ = cond_prob_pi_matrix(X, pi, 0);
    [sZet, protr] = seed_clusters(piZ, protr, min_class, 0);
    [rA,~,~] = fisher_set(X, rZet,alpha, 0);
    if isequal(Zloglik,rZet) || (r == maxR) %  
      break
    else
      Zloglik = rZet;
    end
  end
  figure, imagesc(rZet), xlabel(t), ylabel(r)
  if isequal(tZet, Zloglik) || (t == maxT)
    break
  else
    tZet = Zloglik;
  end
end
% delete(gcp('nocreate'))
figure, imagesc(tZet)
