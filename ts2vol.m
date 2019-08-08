function bold=ts2vol(bold2D,mask)
% bold2D must be nTR x (nx x ny x nz)
% TIME DIMENSION MUST BE FIRST
% bold is nx x ny x nz x nTR
if nargin<2, error('ts2vol needs two arguments'); end
mask=logical(mask);
[nx,ny,nz]=size(mask);
nTR=size(bold2D,1);

bold=zeros(nTR,nx*ny*nz);
bold(:,mask)=bold2D;
bold=reshape(bold,[nTR nx ny nz]);
bold=permute(bold,[2 3 4 1]);
