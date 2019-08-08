clear all; close all;

basePath='../data/S06/NII';
boldFilename=fullfile(basePath,'oBoldEcho1+tlrc.BRIK');
brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');
roiMaskFilename=fullfile(basePath,'resampled_roi+tlrc.BRIK');

[~, bold, iInfo, ~] = BrikLoad (boldFilename);
[~, brainMask, iInfo, ~] = BrikLoad (brainMaskFilename);
[~, roiMask, iInfo, ~] = BrikLoad (roiMaskFilename);
%%
% look at single slice
figure;
subplot(2,2,1);
imagesc(bold(:,:,40,1)); axis on; colormap bone
title('BOLD signal at a given slice for time sample 1');

% look at time series for a given voxel
ts=squeeze(bold(32,38,30,:));
subplot(2,2,2); plot(ts);
title('BOLD time series for a single voxel');

% mask out the brain
subplot(2,2,3);
imagesc(brainMask(:,:,40)); colormap bone
title('Brain Mask');

% go from volume to 2D time series
brainBoldTimeSeries=vol2ts(bold,logical(brainMask));
brainBold4D=ts2vol(brainBoldTimeSeries,logical(brainMask));

%% plot the first 100 voxels in the brain
figure
plot(brainBoldTimeSeries(:,1:10));
legend

%% extract the data in the ROI and in the brain mask
roiMask=logical(roiMask);
brainMask=logical(brainMask);
finalMask=roiMask>0 & brainMask>0;
roiBold=vol2ts(bold,finalMask);

%%
figure;
plot(roiBold(:,[1 250]))
% correlate two voxels with one another
rho=corr(roiBold);
figure; imagesc(rho);




