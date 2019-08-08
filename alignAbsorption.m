% 08/28/18: this doesn't work (ABANDON)
% averaging the absorption is not appropriate due to very small region of
% support of non-negligible absorption.  the new strategy is to project the
% absorption onto the average head.  this will require defining the origin
% of the laser on the standard head
%
% take the A-date generated in preprocessBold and average them, while also
% aligning them to match what was done to the BOLD
clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

% constants
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);
Aprefix=['A-' date];
radMask=19.2; depthMask=26.5; % original values
roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
muRoiPrefix=['mu_' roiPrefix];
nX=64; nY=76; nZ=60; % dimensions of BOLD

%%
% first pass through data to get the mean centroid
centroid=zeros(3,nSubjects);
for s=1:nSubjects
    s
    subjIndx=subjStrs{s};
    basePath=['../data/' subjIndx '/NII/'];
    brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');
    roiMaskFilename=fullfile(basePath,['resampled_' roiPrefix '+tlrc.BRIK']);
    muRoiMaskFilename=fullfile(basePath,['resampled_' muRoiPrefix '+tlrc.BRIK']);
    [~, roiMask, ~, ~] = BrikLoad (roiMaskFilename);
    if isempty(roiMask)
        [~, roiMask, ~, ~] = BrikLoad (muRoiMaskFilename);
    end
    roiMask=logical(roiMask);
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    finalMask=roiMask>0 & brainMask>0;
    [X,Y,Z]=ndgrid(1:nX,1:nY,1:nZ);
    xmask=X(finalMask); ymask=Y(finalMask); zmask=Z(finalMask);
    centroid(:,s)=[mean(xmask);mean(ymask);mean(zmask)];
end
muCentroid=mean(centroid,2); % try to align everyone's bold to for a coherent spatial average
cshifts=round(centroid-repmat(muCentroid,1,nSubjects));

%%
% align the absorption the same way that the BOLD was
muA=zeros(nX,nY,nZ);
for s=3:numel(subjStrs) % start at 3 because 1-2 don't have absorption
    s
    subjIndx=subjStrs{s};
    basePath=['../data/' subjIndx '/NII/'];
    absFilename=fullfile(basePath,['resampled_' Aprefix '+tlrc.BRIK']);
    [~, thisA, aInfo, ~] = BrikLoad (absFilename);
    thisAs=circshift(thisA,[cshifts(1,s) cshifts(2,s) cshifts(3,s)]);
    muA=muA+thisAs/nSubjects;
end

%%

% write out the average bold here
%[~, inputBold, info, ~] = BrikLoad (inputBoldFilename); % just to get "Info"
aInfo.TypeName='float';
aInfo.TypeBytes=4;
aInfo.BRICK_TYPES=3;
[err,ErrMessage,info]=WriteBrikWrap('../data/SAVG/NII/',muA,aInfo,['resampled_mu' Aprefix],'tlrc');

