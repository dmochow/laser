clear all; close all; clc
% construct a mask of the laser absorption area to pass as an
% "anat_follower' to afni_proc.py

subjStr='S23';
pathToData=['../data/' subjStr];
origPath=pwd;
laserOriginPrefix='laserOriginV';
anatFilename=fullfile(pathToData,'anat+orig');
anatWithSkullFilename=fullfile(pathToData,'anatWithSkull+orig');
brainMaskFilename=fullfile(pathToData,'resampled_brain_mask+tlrc.BRIK'); % note that we are working with functional data in the talairach space here
resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
anatBrainMaskFilename=fullfile(pathToData,'brain_mask+orig'); 
biopacFilename=fullfile(pathToData,'subj0024_032819b.mat');
radMask=21.3; depthMask=39.2; % 99 percent min area


% ROI GENERATION
[urf,ulf,urb,lrf,midsag]=getMarkerCoords(subjStr);
laserOrigin = getLaserOrigin(ulf,urf,urb);
laserOrigin = fineTuneOrigin(anatWithSkullFilename,laserOrigin,ulf,urf,lrf);
save(fullfile(pathToData,'laserOrigin.mat'),'laserOrigin');
roiMask = projectLaser(anatBrainMaskFilename,laserOrigin,urf,ulf,urb,lrf,radMask,depthMask);  % project the real laser into the mri
[~, ~, Info, ~] = BrikLoad (anatFilename); % only to get the 'Info' field
roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
WriteBrikWrap(pathToData,roiMask,Info,roiPrefix,'orig');