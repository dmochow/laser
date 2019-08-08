% look at individual time courses
% compare the means and standard deviations across the 3 segments
clear all; close all; clc

!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

% the three time segments
t1=[1:215];
t2=[216:430];
t3=[431:645];

echoStr='1';
filename=['shs8oBoldEcho' echoStr '+tlrc.BRIK'];

for s=1:nSubjects
    pathToData=['../data/' subjStrs{s} '/NII/']
    
    [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,filename));
    
    brainMaskFilename='resampled_brain_mask+tlrc';
    [~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
    brainMask=logical(brainMask);
    
    roiMaskFilename='resampled_roi_r21_z39+tlrc.BRIK';
    [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
    roiMask=logical(roiMask);
    
    if isempty(roiMask) % not all subjects had markers
        roiMaskFilename='resampled_mu_roi_r21_z39+tlrc.BRIK';
        [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
        roiMask=logical(roiMask);
    end
    
    finalMask=roiMask&brainMask;
    
    %%
    ts=vol2ts(bold,finalMask);
    muTs=mean(ts,2); % mean across all voxels in roi
    
    means(s,:)=[mean(muTs(t1)) mean(muTs(t2)) mean(muTs(t3))];
    powers(s,:)=[std(muTs(t1)) std(muTs(t2)) std(muTs(t3))];
    
end
%%
figure; 
subplot(221); 
bar([1 2 3],mean(means));
subplot(222);
bar([1 2 3],mean(powers));



