clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

% make a mask for all subjects without markers by averaging all other ROIs

% all subjects with a mask
subjStr={'S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19'};
nSubjects=numel(subjStr);

maskFilename='resampled_control_roi+tlrc';
allMasks=cell(nSubjects,1);
nVoxels=zeros(nSubjects,1);

for s=1:nSubjects
    s
    path=['../data/' subjStr{s} '/NII/'];
    [err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,maskFilename));
    
    allMasks{s}=mask;
    nVoxels(s)=sum(mask(:)>0);
    
end

%%
muMask=mean(cat(4,allMasks{:}),4);
muNvoxels=round(mean(nVoxels));
[muMaskVals,sortind]=sort(muMask(:),'descend');
keepInds=sortind(1:muNvoxels);

muMaskFinal=zeros(size(muMask));
muMaskFinal(keepInds)=1;

% find the set of voxels such that muNvoxels
%%
% figure;
% for i=1:60
% imagesc(muMaskFinal(:,:,i)); colormap bone
% title(i);
% drawnow;
% pause;
% end

%%
%% write out average mask to all subjects
fullSubjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','SAVG'};
nFullSubjects=numel(fullSubjStr);
for s=1:nFullSubjects
    s
    path=['../data/' fullSubjStr{s} '/NII/'];
    
    InfoWrite=Info;
    opt.Scale=0;
    opt.Prefix='mucontrolroi';
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (muMaskFinal, InfoWrite,opt);
    cd(origPath);
    
end



