clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

% project from laser origin along LOS to brain

!export PATH=$PATH:/Users/jacekdmochowski/abin
% params go here
subjStr='S10';
ulf=[141 23 181]+[1 1 1]; 
urf=[158 28 177]+[1 1 1]; 
urb=[165 9 180]+[1 1 1]; 
lrf=[155 18 137]+[1 1 1];

dLeft=15;
dRight=15;

% subjStr='S11';
% ulf=[127 48 197]+[1 1 1]; 
% urf=[148 55 193]+[1 1 1]; 
% urb=[161 25 213]+[1 1 1]; 
% lrf=[149 30 158]+[1 1 1];
% NB
% dim 1: sagittal
% dim 2: coronal
% dim 3: axial

% try to compute brain mask in afni
% THIS DOESN'T WORK DUE TO MATLAB NOT SEEING THE BASH PATH
% pathToData=['../data/' subjStr '/NII/'];
% origPath=pwd;
% cd(pathToData);
% str='!3dSkullStrip anat.nii';
% eval(str);
% cd(origPath);

% processing begins here
niiFilename=['../data/' subjStr '/NII/anat.nii'];
brainMaskFilename=['../data/' subjStr '/NII/brain_mask+orig'];
outNiifilename=['../data/' subjStr '/NII/roiMask.nii'];

% project from origin, along direction vector, until we are in the brain
% mask
[err, mask, Info, ErrMessage] = BrikLoad (brainMaskFilename);
laserOrigin = getLaserOrigin(ulf,urf,urb); % origin in afni space!

%%
mridim=size(mask);

%% 
% run the search
dr=0.1; % radial increment
currentPoint=laserOrigin;
maxIter=1e6;

crossprod=cross((ulf-urf),(ulf-lrf));
crossprod=crossprod/norm(crossprod);

% if pointing to the front (anterior), flip it
if crossprod(2)<0
    crossprod(2)=-crossprod(2);
end

for i=1:maxIter
    currentPoint=currentPoint+dr*crossprod;
    xq=round(currentPoint); % query point which must be discrete!
    if mask(xq(1)+1,mridim(2)-xq(2),xq(3)+1)
        break
    else
        continue
    end
end

nearestBrainPoint=xq;

%% now figure out the parallel planes on each side of the nearestBrainPoint
lrvec=urf-ulf;
lrvec=lrvec/norm(lrvec);
if lrvec(1)>0
    leftNearestBrainPoint=nearestBrainPoint-dLeft*lrvec; % move left 
    rightNearestBrainPoint=nearestBrainPoint+dLeft*lrvec; 
else
    leftNearestBrainPoint=nearestBrainPoint+dLeft*lrvec; 
    rightNearestBrainPoint=nearestBrainPoint-dLeft*lrvec; 
end

nearestBrainPoint
leftNearestBrainPoint
rightNearestBrainPoint


%% parallel planes on top/bottom of nearestBrainPoint

%%
% for debugging the coordinate transformation between Afni and Matlab via
% BrikLoad
% figure; 
% subplot(221);
% imagesc(squeeze(mask(x+1,:,:)));
% subplot(222);
% imagesc(squeeze(mask(:,mridim(2)-y,:)));  % don't ask, it works
% subplot(223);
% imagesc(squeeze(mask(:,:,z+1)));





