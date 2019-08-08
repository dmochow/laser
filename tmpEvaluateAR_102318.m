

myPath=['../data/SAVG/NII'];
thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
thisBefAftFilename=['thresh_isSigBefAft_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
[~,isSigBefDur,~,~]=BrikLoad(fullfile(myPath,thisBefDurFilename));
[~,isSigBefAft,~,~]=BrikLoad(fullfile(myPath,thisBefAftFilename));
fs=1/2.8;
coeffsPreStim=arIndCoeffs{end}(:,:,[1 2],2);
coeffsPrePost=arIndCoeffs{end}(:,:,[1 3],2);

nfft=32;

freqs=(0:nfft-1)/nfft*fs;
isSig1DbefDur=vol2ts(isSigBefDur,brainMask);

coeffsPreStim=-coeffsPreStim; %
coeffsPreStim=cat(1,ones(1,size(coeffsPreStim,2),size(coeffsPreStim,3)),coeffsPreStim);
A=fft(coeffsPreStim,nfft,1);
C=1./(A+eps); %invert to get the AR spectrum
absC=abs(C);
tmpSbefDur=squeeze(mean(absC(:,isSig1DbefDur==1,:),2));
tmpNbefDur=squeeze(mean(absC(:,isSig1DbefDur==0,:),2));

% before vs after
coeffsPrePost=-coeffsPrePost; %
coeffsPrePost=cat(1,ones(1,size(coeffsPrePost,2),size(coeffsPrePost,3)),coeffsPrePost);
A=fft(coeffsPrePost,nfft,1);
C=1./(A+eps); %invert to get the AR spectrum
absC=abs(C);
%     tmpSbefAft=squeeze(mean(absC(:,isSig1DbefAft==1,:),2));
%     tmpNbefAft=squeeze(mean(absC(:,isSig1DbefAft==0,:),2));

% 09/21/18: substitute which voxels we shoe
tmpSbefAft=squeeze(mean(absC(:,isSig1DbefDur==1,:),2));
tmpNbefAft=squeeze(mean(absC(:,isSig1DbefDur==0,:),2));


%%
hs(1)=subplot(1,3,1);
hstem=stem(freqs,[tmpSbefDur(:,1) tmpSbefDur(:,2) tmpSbefAft(:,2)],'filled');
set(hstem(3),'MarkerEdgeColor',[0.7 0.7 0.7]);
set(hstem(3),'MarkerFaceColor',[0.7 0.7 0.7]);
%axis square
xlim([0 fs/2])