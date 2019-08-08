clear all; close all; clc
addpath(genpath('../../COMMON'));

%
subjStr='S10';

niiFilename=['../data/' subjStr '/NII/anat.nii'];
outNiifilename=['../data/' subjStr '/NII/roiMask.nii'];

% NB
% dim 1: sagittal
% dim 2: coronal
% dim 3: axial

%% get front face
% % S10
ulf=[141 23 181]+[1 1 1]; 
urf=[158 28 177]+[1 1 1]; 
urb=[165 9 180]+[1 1 1]; 

% S11
% ulf=[127 48 197]+[1 1 1]; 
% urf=[148 55 193]+[1 1 1]; 
% urb=[161 25 213]+[1 1 1]; 

% define front face
crossprod=cross((ulf-urf),(ulf-urb));
crossprod=crossprod/norm(crossprod);

% define left-right vector
lrvec=urf-ulf;
lrvec=lrvec/norm(lrvec);

dDown=20; % voxels down from urf
dLeft=10; % voxels left from urf-dDown

if crossprod(3)>0
    cc=urf-dDown*crossprod; % move down
else
    cc=urf+dDown*crossprod; 
end

if lrvec(1)>0
    cc=cc-dLeft*lrvec; % move left 
else
    cc=cc+dLeft*lrvec; 
end

cc


% %%
% nii=load_untouch_nii(niiFilename);
% img=nii.img;
% nx=size(img,1);
% ny=size(img,2);
% nz=size(img,3);
% [X,Y,Z]=ndgrid(1:nx,1:ny,1:nz);
% %x1(2)=256-x1(2); x2(2)=256-x2(2); x3(2)=256-x3(2);


