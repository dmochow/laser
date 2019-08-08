clear all; close all; clc
addpath(genpath('../../COMMON'));

%
subjStr='S07';

niiFilename=['../data/' subjStr '/NII/anat.nii'];
outNiifilename=['../data/' subjStr '/NII/roiMask.nii'];

% NB
% dim 1: sagittal
% dim 2: coronal
% dim 3: axial

% internal parameters
beamHeight=20;
llf=[125 17 162]+[1 1 1]; % adding 1s because AFNI indexes from 0
ulf=[127 16 171]+[1 1 1];
ulb=[131 2 175]+[1 1 1];
lrf=[156 26 157]+[1 1 1]; 
urf=[158 27 166]+[1 1 1];
urb=[165 7 168]+[1 1 1];
isRightPolarity=1;
isLeftPolarity=1;
isBelowPolarity=1;
isAbovePolarity=0;

rangeCutoff=50;  % how close to center of face to allow ROI



%
nii=load_untouch_nii(niiFilename);
img=nii.img;

%
nx=size(img,1);
ny=size(img,2);
nz=size(img,3);
[X,Y,Z]=ndgrid(1:nx,1:ny,1:nz);

%
isRight = makeRoiMask(ulf,llf,ulb,X,Y,Z,isRightPolarity);
isLeft = makeRoiMask(lrf,urf,urb,X,Y,Z,isLeftPolarity);

%% define "below" mask
x1=ulf;
x2=ulb;
x3=urb;
y1=x1; y1(2)=256-y1(2);
y2=x2; y2(2)=256-y2(2);
y3=x3; y3(2)=256-y3(2);
cross_u=cross(y1-y2,y1-y3);
cross_u=cross_u/norm(cross_u);
if cross_u(3)>0 % pointing up
    deld=beamHeight;
else
    deld=-beamHeight;     % pointing down
end
z1=y1+deld*cross_u;
z2=y2+deld*cross_u;
z3=y3+deld*cross_u;
q1=z1; q1(2)=256-z1(2);
q2=z2; q2(2)=256-z2(2);
q3=z3; q3(2)=256-z3(2);
isBelow = makeRoiMask(q1,q2,q3,X,Y,Z,isBelowPolarity);

%% define "above" mask
% here derived as parallel to the below mask
x1=ulf;
x2=ulb;
x3=urb;
y1=x1; y1(2)=256-y1(2);
y2=x2; y2(2)=256-y2(2);
y3=x3; y3(2)=256-y3(2);
cross_u=cross(y1-y2,y1-y3);
cross_u=cross_u/norm(cross_u);
if cross_u(3)>0 % pointing up
    deld=-beamHeight;
else
    deld=beamHeight;     % pointing down
end
z1=y1+deld*cross_u;
z2=y2+deld*cross_u;
z3=y3+deld*cross_u;
q1=z1; q1(2)=256-z1(2);
q2=z2; q2(2)=256-z2(2);
q3=z3; q3(2)=256-z3(2);
isAbove = makeRoiMask(q1,q2,q3,X,Y,Z,0);

%% range mask
faceCenter=(llf+ulf+lrf+urf)/4;
faceCenter(2)=256-faceCenter(2);
FCx=repmat(faceCenter(1),[size(X,1) size(X,2) size(X,3)]);
FCy=repmat(faceCenter(2),[size(Y,1) size(Y,2) size(Y,3)]);
FCz=repmat(faceCenter(3),[size(Z,1) size(Z,2) size(Z,3)]);
RANGE= sqrt( (X-FCx).^2 + (Y-FCy).^2 + (Z-FCz).^2 );
maskRange=RANGE<rangeCutoff;

%% make mask
mask=zeros(size(img));
mask(and(and(and(isLeft,isRight),isAbove),isBelow))=255;
mask=and(mask,maskRange);

% write mask to nifti
nii2=nii;
nii2.img=mask;
save_untouch_nii(nii2,outNiifilename);


%%
% saving marker coordinates for SDB_ here
% llf=[130 10 163]+[1 1 1]; % adding 1s because AFNI indexes from 0
% ulf=[130 11 171]+[1 1 1];
% ulb=[133 1 174]+[1 1 1];
% lrf=[161 20 161]+[1 1 1]; 
% urf=[161 22 168]+[1 1 1];
% urb=[167 1 172]+[1 1 1];

