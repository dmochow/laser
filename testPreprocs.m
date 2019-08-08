%% 05.21.19
% trying to find the magic preprocessing chain yet again
% further burrowing down the rabbit hole
% in search of that carrot
clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
subjStr='S09';
TR=2.8;
nTRs=645; % hard-coded but needs to be changed if this script is to be adapted for the toggling study
alpha=0.05; % significance level
gmIndex=2; % this will never change
roiMaskFilename='resampled_roi_r21_z39+tlrc.BRIK';
altRoiMaskFilename='resampled_mu_roi_r21_z39+tlrc.BRIK'; % no markers
segFilename='resampled_Classes+tlrc.BRIK';
boldPrefix='smoothedOutputBoldEcho'; % use this to test different preprocs
biopacFilename='biopac.mat';

pathToData=['../data/myoutput/' subjStr '/'];
pathToBiopac=['../data/scanner_data/' subjStr '/'];

% illumination ROI
[~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
roiMask=logical(roiMask);

% grey matter mask
resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
[~, seg, ~, ~] = BrikLoad (resampledSegFilename);
greyMask=logical(seg==gmIndex);

finalMask=roiMask&greyMask;

% TODO: load biopac here
% define onset and offset time
load(fullfile(pathToBiopac,biopacFilename),'data');
[~,maxindsample]=max(diff(data(:,1)));
onsetTimeSec=maxindsample/1000;
onsetTimeTR=round(onsetTimeSec/TR);
[~,minindsample]=min(diff(data(:,1)));
offsetTimeSec=minindsample/1000;
offsetTimeTR=round(offsetTimeSec/TR);

t2=[onsetTimeTR:offsetTimeTR-1];
t3=[offsetTimeTR:nTRs];
X=zeros(nTRs,2);
X(t2,1)=1;
X(t3,2)=1;
Xnull=ones(nTRs,1);

for e=1:3
    echoStr=num2str(e);
    %filename=[boldPrefix echoStr '+tlrc.BRIK'];
    filename=['dsbold_e' echoStr '_tlrc_al+tlrc.BRIK'];
    [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,filename));
    
    %%
    ts=vol2ts(bold,finalMask);
    %ts2=bsxfun(@minus,ts,mean(ts));
    %ts3=bsxfun(@rdivide,ts2,mean(ts));
    
    
    % regress out motion
    motionFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_vr_motion.1D']);
    fid=fopen(motionFilename);
    mals=textscan(fid,'%f %f %f %f %f %f');
    mals=cell2mat(mals);
    dmals=cat(1,zeros(1,6),diff(mals));
    mals=cat(2,mals,dmals);
    mals=mals-repmat(mean(mals),size(mals,1),1);
    ts4 = regressOut(ts.',mals.',1).';
 
   
    
end

figure; 
plot(mean(ts4,2));

