% 10.16.18
% show the AR filtered BOLD and compute variance explained
% percent signal change

clear all; close all; clc
PRINT_FIGURE=0;
nTRs=645;
LASERONSETTR=215;
LASEROFFSETTR=430;
testInds=setdiff(1:nTRs,[1:LASERONSETTR  LASEROFFSETTR:nTRs]);
nPerms=1000;
nTRsToSkip=5; % while gradients stabilize
TR=2.8;
time=(0:nTRs-1)*TR;
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
muBoldFilenames={'s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    's8muwBoldEcho2-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    's8muwBoldEcho3-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'};
filenames={'thresh_isSigBefDur_NLAGS5_s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    'thresh_isSigBefDur_NLAGS5_s8muwBoldEcho2-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    'thresh_isSigBefDur_NLAGS5_s8muwBoldEcho3-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'};
arFilenames={'arCoeffsResidsEcho1-NLAGS5_CONSTANT0.mat','arCoeffsResidsEcho2-NLAGS5_CONSTANT0.mat','arCoeffsResidsEcho3-NLAGS5_CONSTANT0.mat'};
coeffIndsDim3=[1;1;1;2]; % 1=pre vs stim; 2=pre vs post
brainMaskFilename='resampled_brain_mask+tlrc.BRIK';
lineColor=[0 0 0];
winsize=20; % time series smoothing parameter
legDelRight=0.075;
legDelUp=0.025;
nFilenames=numel(filenames);
laserOrigin = [40 1 34]+[1 1 1]; % LR, AP, IS
nfft=32; % for AR spectra
fs=1/TR; % for AR spectra
colors=[0 0.4470 0.7410;0.8500    0.3250    0.0980];

[~,brainMask,~,~]=BrikLoad(fullfile(pathToData,brainMaskFilename));
brainMask=logical(brainMask);

hf=figure;
for e=1:nFilenames
    [~,muBold,~,~]=BrikLoad(fullfile(pathToData,muBoldFilenames{e}));
    muBold(:,:,:,1:nTRsToSkip)=0;
    muBoldTs=vol2ts(muBold,brainMask);
    
    thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    thisBefAftFilename=['thresh_isSigBefAft_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    
    [~,isSigBefDur,~,~]=BrikLoad(fullfile(pathToData,thisBefDurFilename));
    %sigBefDurTs=vol2ts(muBold,isSigBefDur>0);
    %muSigBefDurTs=mean(sigBefDurTs,2);
    
    
    % AR coefficients
    load(fullfile(pathToData,arFilenames{e}),'coeffsPreStim','coeffsPrePost');
    
    % before vs during
    coeffsPreStim=-coeffsPreStim; % 
    coeffsPreStim=cat(1,ones(1,size(coeffsPreStim,2),size(coeffsPreStim,3)),coeffsPreStim);
    
    %%
    for v=1:size(muBoldTs,2)
        v
        fBold1(:,v)=filter(1,coeffsPreStim(:,v,1),muBoldTs(:,v));
        fBold2(:,v)=filter(1,coeffsPreStim(:,v,2),muBoldTs(:,v));
    end
    fBold1_v=ts2vol(fBold1,brainMask);
    fBold2_v=ts2vol(fBold2,brainMask);
    sigfBold1Ts=vol2ts(fBold1_v,isSigBefDur>0);
    mu1=mean(sigfBold1Ts,2);
    
    %%
    
    A=fft(coeffsPreStim,nfft,1);
    C=1./(A+eps); %invert to get the AR spectrum
    absC=abs(C);

    % before vs after
    coeffsPrePost=-coeffsPrePost; % 
    coeffsPrePost=cat(1,ones(1,size(coeffsPrePost,2),size(coeffsPrePost,3)),coeffsPrePost);
    A=fft(coeffsPrePost,nfft,1);
    C=1./(A+eps); %invert to get the AR spectrum
    absC=abs(C);
  
end

