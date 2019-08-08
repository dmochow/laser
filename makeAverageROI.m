% make a mask for all subjects without markers by averaging all other ROIs
%
% 07/25/18  account for roi suffixes and include generation of averaged
% control ROI, making makeAverageControlROI.m obsolete
% average in anatomical talairach space, then resample on bold grid
%
% 05/03/19 adapting to new analysis pipeline (i.e., runAfniProc)
% this script should probably be moved to runAfniProc somehow
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

% here we specify the ROI prefix
%radMask=19.2; depthMask=26.5; % original values
radMask=21.3; depthMask=39.2; % min area 99
%roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
roiPrefix='follow_ROI_LASER';
roiFilename=[roiPrefix '+tlrc'];
muRoiPrefix=['mu_' roiPrefix];

% all subjects with markers to indicate laser position
subjStr={'S04','S05','S06','S07','S09','S10','S11','S12','S13','S14',...
    'S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStr);
allMasks=cell(nSubjects,1); allControlMasks=cell(nSubjects,1);
nVoxels=zeros(nSubjects,1); nControlVoxels=zeros(nSubjects,1);

for s=1:nSubjects
    s
    path=['../data/output/' subjStr{s} '/' subjStr{s} '.results/'];
    % true laser
    [err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,roiFilename));
    allMasks{s}=mask;
    nVoxels(s)=sum(mask(:)>0);
end

%%
% true laser
muMask=mean(cat(4,allMasks{:}),4);
muNvoxels=round(mean(nVoxels));
[muMaskVals,sortind]=sort(muMask(:),'descend');
keepInds=sortind(1:muNvoxels);
muMaskFinal=zeros(size(muMask));
muMaskFinal(keepInds)=1;

% write out average mask to all subjects without markers
fullSubjStr={'S02','S03'};
nFullSubjects=numel(fullSubjStr);
for s=1:nFullSubjects
    
    path=['../data/output/' fullSubjStr{s} '/' fullSubjStr{s} '.results/'];
    mkdir(path);
    
    % if the files exist, delete them first
    filenames=dir(path);
    for f=1:numel(filenames)
        if (strcmp(filenames(f).name,[muRoiPrefix '+tlrc.BRIK']) || ...
                strcmp(filenames(f).name,[muRoiPrefix '+tlrc.HEAD']) || ...
                strcmp(filenames(f).name,['resampled_' muRoiPrefix '+tlrc.BRIK']) || ...
                strcmp(filenames(f).name,['resampled_' muRoiPrefix '+tlrc.HEAD'])  )
            eval(['delete ' fullfile(path,filenames(f).name)]);
        end
    end
    
    InfoWrite=Info;
    opt.Scale=0;
    opt.Prefix=muRoiPrefix;
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (muMaskFinal, InfoWrite,opt);
    cd(origPath);
    
    % now resample, if a master is available
    cd(path);
    %bashStr3=['!3dresample -master all_runs.' subjStr '+tlrc -prefix resampled_Classes -input Classes+tlrc']; eval(bashStr3);
    if exist(['all_runs.' subjStr{s} '+tlrc.BRIK'])
        str1=['!3dresample -master all_runs.' subjStr{s} '+tlrc -prefix resampled_' muRoiPrefix ' -input ' muRoiPrefix '+tlrc']; eval(str1);  % resample on bold grid
    end
    cd(origPath);
    
end



