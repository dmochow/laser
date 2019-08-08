function stats=myGLM(Y,X,addOnes)
% Y: 1-D BOLD time series
% X: design matrix where number of rows matches that of Y (columns are
% predictors)
% addOnes (defaults to 1): whether to add a column of ones to design matrix
if nargin<3, addOnes=1; end
if nargin<2, error('Need at least two arguments'); end
nTRs=size(Y,1);
if addOnes, X=cat(2,X,ones(nTRs,1)); end

B=pinv(X)*Y;
sse=sum((Y-X*B).^2);
dof=nTRs-size(X,2);

stats.B=B;
stats.sse=sse;
stats.dof=dof;