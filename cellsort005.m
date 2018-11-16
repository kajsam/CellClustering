function cellsort005(X, Zclust, sAr, ticks, ticklabels)

% Sort by cell class
merg = [1 2];

Xclass = X(logical(sum(Zclust(:,merg),2)),:);
[nclass,d] = size(Xclass);
Xclass(1:sum(Zclust(:,merg(1))),:) = X(Zclust(:,merg(1)),:);
Xclass(sum(Zclust(:,merg(1)))+1:nclass,:) = X(Zclust(:,merg(2)),:);
X(logical(sum(Zclust(:,merg),2)),:) = Xclass;

Aclass = sAr(logical(sum(Zclust(:,merg),2)),:);
Aclass(1:sum(Zclust(:,merg(1))),:) = sAr(Zclust(:,merg(1)),:);
Aclass(sum(Zclust(:,merg(1)))+1:nclass,:) = sAr(Zclust(:,merg(2)),:);
sAr(logical(sum(Zclust(:,merg),2)),:) = Aclass;
figure(29), imagesc(sAr)

Zclass = false(nclass,2);
Zclass(1:sum(Zclust(:,merg(1))),merg(2)) = true;
Zclass(sum(Zclust(:,merg(1)))+1:nclass,merg(1)) = true;

cZclust = Zclust;
cZclust(logical(sum(Zclust(:,merg),2)),:) = Zclass;

col = [0.25 0.4 0.7 0.9];
rep = 100;
imZclust = [col(1)*repmat(cZclust(:,1), 1,rep) col(2)*repmat(cZclust(:,2), 1,rep) ...
            col(3)*repmat(cZclust(:,3), 1,rep) col(4)*repmat(cZclust(:,4), 1,rep)];

        
rowA = unique(sAr,'rows','stable'); % Finding replicates. 
rowA = rowA([1 3 2 4],:);
figure(27), imagesc(rowA), colormap(gray), title('Fisher unique')
      
imA = [col(4)*repmat(rowA(1,:),25,1); col(3)*repmat(rowA(2,:),25,1);...
       col(2)*repmat(rowA(3,:),25,1); col(1)*repmat(rowA(4,:),25,1)];
   
imA = [imA zeros(100,4*rep)];
imA = ones(size(imA))*0.1 + imA;

imZclust = ones(size(imZclust))*0.1 + imZclust;

figure(32), imagesc(imA)
figure(30), imagesc([imA; 0.6*X imZclust]), colormap(jet(64))
xticks(ticks)
xticklabels(ticklabels)
xtickangle(90)
title('P3CL 5% filtering')
xlabel('Fisher hypothesis tests, p-value = 0.05, 10637 rejected H0`s (98.6%)')




