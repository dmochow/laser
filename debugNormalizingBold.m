clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%% PREAMBLE
subjStr='S27';
pathToData=['../data/' subjStr];
biopacFilename=fullfile(pathToData,'subj0027_040419.mat');
origPath=pwd;
laserOriginPrefix='laserOriginV';
anatFilename=fullfile(pathToData,'anat+orig');
anatWithSkullFilename=fullfile(pathToData,'anatWithSkull+orig');
brainMaskFilename=fullfile(pathToData,'resampled_brain_mask+tlrc.BRIK'); % note that we are working with functional data in the talairach space here
resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
anatBrainMaskFilename=fullfile(pathToData,'brain_mask+orig');
TR=2.8;
fs=1/TR;
nX=64; nY=76; nZ=60; % dimensions of BOLD
nTR=645; % number of TRs
dummyMask=true(nX,nY,nZ); % used for going from 2D to 4D
Kwm=3; % number of WM PCs to remove
fl=0.1; % low-pass frequency cutoff
nPad=30; % number of samples to pad to avoid filter transient (don't change)
[bf,af]=butter(3,fl/(fs/2)); % compute the filter coefficients here
WM_INDEX=3; % segmentation label for white matter
radMask=21.3; depthMask=39.2; % 99 percent min area

% load in brain mask
[~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
brainMask=logical(brainMask);

%%
% load input bold ts
for e=1:3
    inputBoldFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_tlrc_al+tlrc.BRIK']);
    [~, inputBold, iInfo, ~] = BrikLoad (inputBoldFilename);
    brainTs=vol2ts(inputBold,brainMask);
    tmp(:,:,e)=brainTs;
    
    % create the motion regressors
    motionFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_vr_motion.1D']);
    fid=fopen(motionFilename);
    mals=textscan(fid,'%f %f %f %f %f %f');
    mals=cell2mat(mals);
    dmals=cat(1,zeros(1,6),diff(mals));
    mals=cat(2,mals,dmals);
    mals=mals-repmat(mean(mals),size(mals,1),1);
    
    % regress out motion
    obrainTs = regressOut(brainTs.',mals.',1).';
    tmp2(:,:,e)=obrainTs;
    outputBold=ts2vol(obrainTs,brainMask);
    
    %     % regress out white matter
    %     [~,seg,~,~] = BrikLoad(resampledSegFilename);
    %     wmMask=seg==WM_INDEX;
    %     wmBoldTs=vol2ts(outputBold,wmMask);
    %     [U,S,V]=svd(wmBoldTs,0);
    %     brainBoldTs=vol2ts(outputBold,brainMask); % only regress in the brain mask
    %     oobrainTs = regressOut(brainBoldTs.',U(:,1:Kwm).',1).';
    %     tmp3(:,:,e)=oobrainTs;
    
    % write out the post-processed bold here
    [err,ErrMessage,Info]=WriteBrikWrap(pathToData,outputBold,iInfo,['unsmoothedOutputBoldEcho' num2str(e)],'tlrc');
    
end

%%
% filtering
% filter
%tmp=cat(1,repmat(obrainBoldTs(1,:),nPad,1),obrainBoldTs);
%tmp4=filter(bf,af,tmp3,[],1);
%obrainBoldTs=tmpOut(nPad+1:end,:);

% % 2D --> 4D
% final4D=ts2vol(obrainBoldTs,brainMask);



%%
% v=30000;
% t=430;
% figure;
% subplot(421);
% plot(squeeze(tmp(t,v,:)));
% subplot(422);
% plot(squeeze(tmp(:,v,:)));
% subplot(423);
% plot(squeeze(tmp2(t,v,:)));
% subplot(424);
% plot(squeeze(tmp2(:,v,:)));
% subplot(425);
% plot(squeeze(tmp3(t,v,:)));
% subplot(426);
% plot(squeeze(tmp3(:,v,:)));
% subplot(427);
% plot(squeeze(tmp4(t,v,:)));
% subplot(428);
% plot(squeeze(tmp4(:,v,:)));



