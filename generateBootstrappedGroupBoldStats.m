%% bootstrapped estimates of group bold

clear all; close all; clc

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);
nx=64; ny=76; nz=60; nTRs=645;
ECHOSTR='1';
boldFilename=['shs8oBoldEcho' ECHOSTR '+tlrc'];
brainMaskFilename='resampled_brain_mask+tlrc.BRIK';
nResamples=100;
nfft=32; posFreqInds=1:nfft/2+1;

%% pull in the significant voxels for this echo
avgPath='../data/SAVG/NII';
thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho' ECHOSTR '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
[~,isSigBefDur,~,~]=BrikLoad(fullfile(avgPath,thisBefDurFilename));

% bring in the brain mask
[~,brainMask,~,~]=BrikLoad(fullfile(avgPath,brainMaskFilename));
brainMask=logical(brainMask);

isSigBefDur1D=vol2ts(isSigBefDur,brainMask);

%%
sigTs=zeros(nTRs,nResamples);
nsigTs=zeros(nTRs,nResamples);
sigMuSpectra=zeros(nfft,3,nResamples); % before/during/after
nsigMuSpectra=zeros(nfft,3,nResamples);

%%

for r=1:nResamples
    tic
    sInds = datasample(1:nSubjects,nSubjects);
    %sInds=1:20;
    groupBold=zeros(nx,ny,nz,nTRs); % initialize/reset
    
    % generate the bootstrapped group bold
    for s=1:nSubjects
        [r,s]
        % group bold time course
        pathToData=['../data/' subjStrs{sInds(s)} '/NII'];
        [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,boldFilename));
        groupBold=groupBold+bold/nSubjects;
    end
    
    % get the time course
    sigGroupBold=vol2ts(groupBold,isSigBefDur);
    nsigGroupBold=vol2ts(groupBold,~isSigBefDur & brainMask);
    sigTs(:,r)=mean(sigGroupBold,2);
    nsigTs(:,r)=mean(nsigGroupBold,2);
    
    % get the AR coefficients and spectrum
    groupBoldTs=vol2ts(groupBold,brainMask);
    [coeffs,arSpectrum,pvals,tstats] = computeArCoeffs(groupBoldTs);
    
    sigSpectra=arSpectrum(:,isSigBefDur1D==1,:);
    nsigSpectra=arSpectrum(:,isSigBefDur1D==0,:);
    
    sigMuSpectra(:,:,r)=squeeze(mean(sigSpectra,2));
    nsigMuSpectra(:,:,r)=squeeze(mean(nsigSpectra,2));
    toc

    save(['../data/precomputed/boostrappedGroupStatsEcho' ECHOSTR '.mat'],'sigTs','nsigTs','sigMuSpectra','nsigMuSpectra');
end

%%

% tmp=vol2ts(muBold,isSigBefDur);
% figure; plot(mean(tmp,2));
