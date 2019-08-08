function [Y,Z]=mynorm(X)

goodind=sum(X.^2)>0.01;
Y=zeros(size(X));
Y(:,goodind)=X(:,goodind)./repmat(sqrt(sum(X(:,goodind).^2)),size(X,1),1);
Z=Y(:,goodind);