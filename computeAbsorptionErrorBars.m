% compute error bars around GM/WM/CSF variability
%
% since we aligned laser position across subjects, assume the same laser
% origin and instead get variability from segmentation differences

clear all; close all; clc

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

%% anatomy
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
brainMaskFilename=fullfile(pathToData,'brain_mask+tlrc.BRIK');
absorptionFilename=fullfile(pathToData,'Amu+tlrc.BRIK');
[~, absorption, ~, ~] = BrikLoad (absorptionFilename);
absorption=absorption(end:-1:1,:,:);

voxelVolume=0.1^3; % cm^3

for s=1:nSubjects
    s
    [~,labels,~,~]=BrikLoad(['../data/' subjStrs{s}  '/NII/Classes+tlrc.BRIK']);
    labels=labels(end:-1:1,:,:);
    
    csfTotal=sum(absorption(labels==1))*voxelVolume;
    gmTotal=sum(absorption(labels==2))*voxelVolume;
    wmTotal=sum(absorption(labels==3))*voxelVolume;
    
    csfMax=max(absorption(labels==1));
    gmMax=max(absorption(labels==2));
    wmMax=max(absorption(labels==3));
    
    csfTotals(s)=csfTotal;
    gmTotals(s)=gmTotal;
    wmTotals(s)=wmTotal;
    
    csfMaxima(s)=csfMax;
    gmMaxima(s)=gmMax;
    wmMaxima(s)=wmMax;
    
end

save('../data/precomputed/absorptionErrorBars','csfTotals','gmTotals','wmTotals','csfMaxima','gmMaxima','wmMaxima');