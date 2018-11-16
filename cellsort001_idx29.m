function [X, Zet, A] = cellsort001_idx29(X, Zet, A)


figure(25), imagesc(Zet)

% Sort by cell class
merg = [3 4];
Xmov = X(logical(sum(Zet(1:250,merg),2)),:);
[nmov,d] = size(Xmov);
Xmov(1:sum(Zet(1:250,merg(1))),:) = X(Zet(1:250,merg(1)),:);
Xmov(sum(Zet(1:250,merg(1)))+1:nmov,:) = X(Zet(1:250,merg(2)),:);
X(logical(sum(Zet(1:250,merg),2)),:) = Xmov;

Amov = A(logical(sum(Zet(1:250,merg),2)),:);
Amov(1:sum(Zet(1:250,merg(1))),:) = A(Zet(1:250,merg(1)),:);
Amov(sum(Zet(1:250,merg(1)))+1:nmov,:) = A(Zet(1:250,merg(2)),:);
A(logical(sum(Zet(1:250,merg),2)),:) = Amov;

Zmov = Zet(logical(sum(Zet(1:250,merg),2)),:);
Zmov(1:sum(Zet(1:250,merg(1))),:) = Zet(Zet(1:250,merg(1)),:);
Zmov(sum(Zet(1:250,merg(1)))+1:nmov,:) = Zet(Zet(1:250,merg(2)),:);
Zet(logical(sum(Zet(1:250,merg),2)),:) = Zmov;
Zet = Zet(:,[1 2 4 3]);

figure(26), imagesc(Zet)





