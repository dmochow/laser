% 10.24.18 compute Fourier transform on raw data
% NB: significant voxels have been taken from mean-bold analysis

clear all; close all; clc

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

% global parameters
nEchos=3;
nTRs=645;
LASERONSETTR=215;
LASEROFFSETTR=430;
nfft=256;
TR=2.8;
fs=1/TR;
freqs=(0:nfft-1)/nfft*fs;

myPath='../data/SAVG/NII';

%%
indSpectra=cell(nSubjects,nEchos); 
for e=1:nEchos
    thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    thisBefAftFilename=['thresh_isSigBefAft_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    [~,isSigBefDur,~,~]=BrikLoad(fullfile(myPath,thisBefDurFilename));
    [~,isSigBefAft,~,~]=BrikLoad(fullfile(myPath,thisBefAftFilename));
    
    for s=1:nSubjects  
        [s,e]
        pathToData=['../data/' subjStrs{s} '/NII'];
        brainMaskFilename='resampled_brain_mask+tlrc';
        [~, brainMask, ~, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
        brainMask=logical(brainMask);
        nVoxels=sum(brainMask(:));

        dataFilename=['s8oBoldEcho' num2str(e) '+tlrc'];
        [~, bold, info, ~] = BrikLoad (fullfile(pathToData,dataFilename));
        
        %boldTs=vol2ts(bold,brainMask);
        boldTs=vol2ts(bold,isSigBefDur);
        
        % spectrogram
        thisSpecGram = spectrogram(mean(boldTs,2),32,24,256);
        allSpecGrams(:,:,s,e)=thisSpecGram;
        
%         boldTs3(:,:,1)=boldTs(1:LASERONSETTR,:);
%         boldTs3(:,:,2)=boldTs(LASERONSETTR+1:LASEROFFSETTR,:);
%         boldTs3(:,:,3)=boldTs(LASEROFFSETTR+1:end,:);
%         
%         indSpectra{s,e}=fft(boldTs3,nfft,1);
%   
%         clear boldTs3;
    end
    
    
end

%%
muSpecGram=squeeze(mean(abs(allSpecGrams),3));
figure
for e=1:nEchos
    subplot(1,3,e);
    imagesc([0 30],[0 0.17],muSpecGram(:,:,e));
end
print -dpng ../figures/rawSpecGrams

%save('../data/precomputed/individualSpectra.mat','indSpectra');
%%
% for s=1:nSubjects
%     for e=1:nEchos
%         thisS=indSpectra{s,e};
%         absS(:,:,e,s)=squeeze( mean ( abs(thisS),2 ) );
%     end
% end
%     
% muAbsS=mean(absS,4);
% %%
% figure
% for e=1:nEchos
%     hs(e)=subplot(1,3,e);
%     plot(freqs,mean(muAbsS(:,:,e),3));
%     xlim([0 fs/2]);
%     axis square;
%     legend('Before','During','After');
% end
% 
% print -dpng ../figures/rawFFTs.png
