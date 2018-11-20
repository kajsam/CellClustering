function [X, Zet, A] = cellsort001_random(X, Zet, A)

% Sort by cell class
% figure, imagesc(Zet), title('Before')

[n,d] = size(X);

merg = [1 2];
mZet = Zet(:,merg);
mZet(1:871,:) = false;
% figure, imagesc(mZet), title('Sort these')
Xmov = X(logical(sum(mZet,2)),:);
[nmov,d] = size(Xmov)
sum(mZet(:,1))
sum(mZet(:,2))
Xmov(1:sum(mZet(:,1)),:) = X(mZet(:,1),:);
Xmov(sum(mZet(:,1))+1:nmov,:) = X(mZet(:,2),:);
X(logical(sum(mZet,2)),:) = Xmov;

Amov = A(logical(sum(mZet,2)),:);
Amov(1:sum(mZet(:,1)),:) = A(mZet(:,1),:);
Amov(sum(mZet(:,1))+1:nmov,:) = A(mZet(:,2),:);
A(logical(sum(mZet,2)),:) = Amov;

Zmov = Zet(logical(sum(mZet,2)),:);
Zmov(1:sum(mZet(:,1)),:) = Zet(mZet(:,1),:);
Zmov(sum(mZet(:,1))+1:nmov,:) = Zet(mZet(:,2),:);
Zet(logical(sum(mZet,2)),:) = Zmov;

Zet = Zet(:,[2 1 4 3]);

% figure, imagesc(Zet), title('After')
