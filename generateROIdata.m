% this script after preprocess bold
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

useAverageROI=0; % if we want to apply the same average ROI to all subjects
subjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19'};
nSubjects=numel(subjStr);
writeBriks=0;

maskFilename='resampled_roi+tlrc';
controlMaskFilename='resampled_control_roi+tlrc';
altMaskFilename='muroi+tlrc'; % if subject doesn't have markers
altControlMaskFilename='muControlRoi+tlrc'; % if subject doesn't have markers
brainMaskFilename='resampled_brain_mask+tlrc';
epiFilename_1='smnudsbold_e1_tlrc_al+tlrc';
epiFilename_2='smnudsbold_e2_tlrc_al+tlrc';
epiFilename_3='smnudsbold_e3_tlrc_al+tlrc';
biopacFilename='biopac.mat';
TEs=[12.8 34.13 55.46];
nEchos=numel(TEs);
TR=2.8;
alpha=0.05;


allOnsets=zeros(nSubjects,2);
allBoldsRoi=cell(nSubjects,1);
allBoldsGlobal=cell(nSubjects,1);
allBoldsControlRoi=cell(nSubjects,1);


%%
for s=1:nSubjects
    s
    path=['../data/' subjStr{s} '/NII/'];
    biopacPath=['../data/' subjStr{s} '/BIOPAC/'];
    
    [err, roiMask, InfoMask, ~] = BrikLoad (fullfile(path,maskFilename));
    if err || useAverageROI
        [~, roiMask, InfoMask, ~] = BrikLoad (fullfile(path,altMaskFilename));
    end
    
    [err, controlRoiMask, InfoMask, ~] = BrikLoad (fullfile(path,controlMaskFilename));
    if err || useAverageROI
        [~, controlRoiMask, InfoMask, ~] = BrikLoad (fullfile(path,altControlMaskFilename));
    end
    
    [~, brainMask, Info, ~] = BrikLoad (fullfile(path,brainMaskFilename));
    [~, bold_1, Info, ~] = BrikLoad (fullfile(path,epiFilename_1));
    [~, bold_2, Info, ~] = BrikLoad (fullfile(path,epiFilename_2));
    [~, bold_3, Info, ~] = BrikLoad (fullfile(path,epiFilename_3));
    
    % compute the mask
    finalMask=roiMask>0 & brainMask>0; % in the laser beam AND in the brain mask
    finalControlMask=controlRoiMask>0 & brainMask>0;
    
    maskInds=find(finalMask);
    controlMaskInds=find(finalControlMask);
    
    % grab the bold in the ROI
    tmp=permute(bold_1,[4 1 2 3]);
    bold(1,:,:)=tmp(:,:);
    
    tmp=permute(bold_2,[4 1 2 3]);
    bold(2,:,:)=tmp(:,:);
    
    tmp=permute(bold_3,[4 1 2 3]);
    bold(3,:,:)=tmp(:,:);
    
    boldMasked=bold(:,:,finalMask);
    allBoldsRoi{s,1}=boldMasked;
    
    boldGlobal=mean(bold,3);
    allBoldsGlobal{s}=boldGlobal;
    
    boldControlMasked=bold(:,:,finalControlMask);
    allBoldsControlRoi{s,1}=boldControlMasked;
    
    % laser onset times from biopac
    load(fullfile(biopacPath,biopacFilename),'data');
    
    [~,maxindsample]=max(diff(data(:,1)));
    onsetTimeSec=maxindsample/1000;
    onsetTimeTR=round(onsetTimeSec/TR);
    
    [~,minindsample]=min(diff(data(:,1)));
    offsetTimeSec=minindsample/1000;
    offsetTimeTR=round(offsetTimeSec/TR);
    
    if offsetTimeTR==onsetTimeTR % error
        offsetTimeTR=onsetTimeTR+214;
    end
    
    allOnsets(s,1)=onsetTimeTR;
    allOnsets(s,2)=offsetTimeTR;
    nTRsLaser=offsetTimeTR-onsetTimeTR+1;
    
end

%%
precomputedFilename=['../data/precomputed/allBoldsRoi ' date];
save(precomputedFilename,'allBoldsRoi','allBoldsGlobal','allOnsets','allBoldsControlRoi','TR','TEs','subjStr','nSubjects','nEchos');



