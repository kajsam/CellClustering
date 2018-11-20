function [X, Zet, A] = cellsort001_idx29(X, Zet, A)

% Sort by cell class
figure, imagesc(Zet), title('Before')

merg = [3 4];
mZet = Zet(:,merg);
mZet(251:end,:) = false;
figure, imagesc(mZet), title('Sort these')
Xmov = X(logical(sum(mZet,2)),:);
[nmov,d] = size(Xmov);
pause
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
Zet = Zet(:,[1 2 4 3]);

figure, imagesc(Zet), title('After')