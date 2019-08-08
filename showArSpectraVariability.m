% show variability estimates of AR spectra across subjects

clear all; close all; clc

% parameters
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);
nEchos=3;

%
load('../data/precomputed/arIndividualData.mat','pvals','tstats','arIndCoeffs');

%pathToData=['../data/' subjStrs{s} '/NII'];


%
for e=1:nEchos
    % grab the significant voxels for this echo
    thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    [~,isSigBefDur,~,~]=BrikLoad(fullfile('../data/SAVG/NII',thisBefDurFilename));
    for s=1:nSubjects
        [s,e]
        pathToData=['../data/' subjStrs{s} '/NII'];
        brainMaskFilename='resampled_brain_mask+tlrc';
        [~, brainMask, ~, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
        thisBrainMask=logical(brainMask);
        theseCoeffs=squeeze(arIndCoeffs{s}(:,:,:,e));  
        for i=1:4
            for j=1:3
                theseTheseCoeffs=squeeze(theseCoeffs(i,:,j));
                theseTheseCoeffsVol=ts2vol(theseTheseCoeffs,thisBrainMask);
                sigCoeffsTs(i,:,j)=vol2ts(theseTheseCoeffsVol,isSigBefDur);
            end
        end
        allSigCoeffsTs{s,e}=sigCoeffsTs;
    end
    clear sigCoeffsTs
end

%%
tmp=cat(4,allSigCoeffsTs{:,1});
tmp=-tmp; %
tmp=cat(1,ones(1,size(tmp,2),size(tmp,3),size(tmp,4)),tmp);
nfft=32;
A=fft(tmp,nfft,1);
C=1./(A+eps); %invert to get the AR spectrum
tmp2=mean(C,4);
tmp3=squeeze(mean(tmp2,2));
%%
figure;
stem(abs(tmp3));
