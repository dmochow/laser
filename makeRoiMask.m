function mask = makeRoiMask(x1,x2,x3,X,Y,Z,gg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%llf=[130 10 163]+[1 1 1]; % adding 1s because AFNI indexes from 0
%ulf=[130 11 171]+[1 1 1];
%ulb=[133 1 174]+[1 1 1];

%[X,Y,Z]=ndgrid(1:nx,1:ny,1:nz);

x1(2)=256-x1(2); x2(2)=256-x2(2); x3(2)=256-x3(2);

% define "left" plane
%cross_l=cross((llf-llb),(llf-ulf));
cross_l=cross((x1-x2),(x1-x3));
xo=x2;
X0=repmat(xo(1),[size(X,1) size(X,2) size(X,3)]);
Y0=repmat(xo(2),[size(Y,1) size(Y,2) size(Y,3)]);
Z0=repmat(xo(3),[size(Z,1) size(Z,2) size(Z,3)]);
if gg
    mask= ( cross_l(1)*(X-X0)+cross_l(2)*(Y-Y0)+cross_l(3)*(Z-Z0) ) > 0;
else
    mask= ( cross_l(1)*(X-X0)+cross_l(2)*(Y-Y0)+cross_l(3)*(Z-Z0) ) < 0;
end
end

