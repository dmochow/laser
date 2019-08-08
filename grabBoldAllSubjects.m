clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

% point to the data
subjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'};
nSubjects=numel(subjStr);
figFilename=['../figures/' num2str(nSubjects) 'subs_grand_average'];

maskFilename='resampled_roi+tlrc';
altMaskFilename='muroi+tlrc';
anatFilename='resampled_anat+tlrc';
epiFilename_1='smnudsbold_e1_tlrc_al+tlrc';
epiFilename_2='smnudsbold_e2_tlrc_al+tlrc';
epiFilename_3='smnudsbold_e3_tlrc_al+tlrc';
% epiFilename_1='dsbold_e1_tlrc_al+tlrc';
% epiFilename_2='dsbold_e2_tlrc_al+tlrc';
% epiFilename_3='dsbold_e3_tlrc_al+tlrc';
biopacFilename='biopac.mat';
TEs=[12.8 34.13 55.46];
nEchos=numel(TEs);

allBolds=cell(nSubjects,1);
allBoldsControl=cell(nSubjects,1);
allMaskInds=cell(nSubjects,1);
allOnsets=zeros(nSubjects,2);

for s=1:nSubjects
    s
    path=['../data/' subjStr{s} '/NII/'];
    biopacPath=['../data/' subjStr{s} '/BIOPAC/'];
    % read in the data
    [err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,maskFilename));
    if err
        [err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,altMaskFilename));
    end
    [err, anat, Info, ErrMessage] = BrikLoad (fullfile(path,anatFilename));
    [err, bold_1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
    [err, bold_2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
    [err, bold_3, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_3));
    
    % compute the mask
    finalMask=mask>0 & anat>0; % in the laser beam AND in the brain
    maskInds=find(finalMask);
    
    % grab the bold in the ROI
    tmp=permute(bold_1,[4 1 2 3]);
    bold(1,:,:)=tmp(:,:);
    
    tmp=permute(bold_2,[4 1 2 3]);
    bold(2,:,:)=tmp(:,:);
    
    tmp=permute(bold_3,[4 1 2 3]);
    bold(3,:,:)=tmp(:,:);
    %
    bold_roi=bold(:,:,finalMask);
    allBolds{s}=bold_roi;
    allMaskInds{s}=maskInds;
    
%     bold_nroi=bold(:,:,~finalMask);
%     allBoldsControl{s}=bold_nroi;
    
    load(fullfile(biopacPath,biopacFilename),'data');
    [~,maxindsample]=max(diff(data(:,1)));
    [~,minindsample]=min(diff(data(:,1)));
    onsetTimeSec=maxindsample/1000;
    onsetTimeTR=round(onsetTimeSec/2.8);
    offsetTimeSec=minindsample/1000;
    offsetTimeTR=round(offsetTimeSec/2.8);
    %[onsetTimeTR,offsetTimeTR]
    allOnsets(s,1)=onsetTimeTR;
    allOnsets(s,2)=offsetTimeTR;
    
    
end

%%
% because I'm eager
time=(0:644)*2.8;
grandBold=mean(cat(3,allBolds{:}),3);
figure; 
subplot(411);
plot(time,grandBold(1,:)); title('Echo 1');
subplot(412);
plot(time,grandBold(2,:)); title('Echo 2');
subplot(413);
plot(time,grandBold(3,:)); title('Echo 3');
subplot(414);
plot(time,mean(grandBold)); title('Echo mean');
print('-depsc','-r600',figFilename);
%%
%save ../data/precomputed/groupRoiData9subjects allBolds allOnsets allMaskInds

% %%
% grandSum=0;
% grandNum=0;
% for s=1:nSubjects
%     tmp=allBoldsControl{s};
%     grandSum=grandSum+sum(tmp,3);
%     grandNum=grandNum+size(tmp,3);
% end
% grandBoldControl=grandSum/grandNum;
% figure; 
% subplot(221);
% plot(grandBoldControl(1,:));
% subplot(222);
% plot(grandBoldControl(2,:));
% subplot(223);
% plot(grandBoldControl(3,:));
% 
% %%
% figure;
% subplot(211);
% plot((0:644)*2.8,mean(grandBold));
% title('ROI');
% subplot(212);
% plot((0:644)*2.8,mean(grandBoldControl));
% title('Outside of ROI');
% 

