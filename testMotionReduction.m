clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%%
TR=2.8;
subjIndx='S11';
basePath=['../data/' subjIndx '/NII/'];
malFile=['../data/' subjIndx '/NII/dsbold_e1_vr_motion.1D'];
inputBoldFilename=[basePath 'dsbold_e1_tlrc_al+tlrc.BRIK'];
brainMaskFilename=[basePath 'resampled_brain_mask+tlrc.BRIK'];
segFilename=[basePath 'resampled_Classes+tlrc.BRIK'];
%roiMaskFilename=[basePath 'muroi+tlrc.BRIK'];
roiMaskFilename=[basePath 'resampled_roi+tlrc.BRIK'];
anatFilename=[basePath 'anat_ns+orig.BRIK'];
motionFilename=[basePath 'dsbold_e1_vr_motion.1D'];

outputBoldFilename=[basePath 'mynudsbold_e1_tlrc_al+tlrc.BRIK'];

[~, inputBold, Info, ~] = BrikLoad (inputBoldFilename);

[~, roiMask, Info, ~] = BrikLoad (roiMaskFilename);
roiMask=logical(roiMask);
[~, brainMask, Info, ~] = BrikLoad (brainMaskFilename);
brainMask=logical(brainMask);
[nx,ny,nz,nTR]=size(inputBold);

%% 3DTPROJECT
str=['delete ' basePath 'mynudsbold_e1_tlrc_al+tlrc.BRIK'];
eval(str);
str=['delete ' basePath 'mynudsbold_e1_tlrc_al+tlrc.HEAD'];
eval(str);
str=['!3dTproject -input ' inputBoldFilename ' -prefix ' basePath 'mynudsbold_e1_tlrc_al -ort ' motionFilename];
eval(str);
[~, outputBold, Info, ~] = BrikLoad (outputBoldFilename);

%% COMPARE WITH REGRESSOUT

iwBoldTs=vol2ts(inputBold);
oBoldTs=vol2ts(outputBold,roiMask);

fid=fopen(motionFilename);
mals=textscan(fid,'%f %f %f %f %f %f');
mals=cell2mat(mals);
mals=mals-repmat(mean(mals),size(mals,1),1);

owBoldTs = regressOut(iwBoldTs.',mals.',1).';
myOutputBold=ts2vol(owBoldTs,logical(ones(nx,ny,nz)));
myoBoldTs=vol2ts(myOutputBold,roiMask);


%myoBoldTs=normc(myoBoldTs);
%[myoBoldTs,Z]=mynorm(myoBoldTs);
% myoBoldTs=myoBoldTs./repmat(sqrt(sum(myoBoldTs.^2,1)),size(myoBoldTs,1),1);


figure;
subplot(211);
plot(mean(oBoldTs,2));
subplot(212);
plot(mean(Z,2));
corrcoef(mean(oBoldTs,2),mean(myoBoldTs,2))


% %% ADD DERIVATIVES AND ABSOLUTE VALUES TO REGRESSION
% iwBoldTs=vol2ts(inputBold);
% oBoldTs=vol2ts(outputBold,roiMask);
% 
% fid=fopen(motionFilename);
% mals=textscan(fid,'%f %f %f %f %f %f');
% mals=cell2mat(mals);
% dmals=cat(1,zeros(1,6),diff(mals));
% amals=abs(mals);
% mals=cat(2,cat(2,mals,dmals),amals);
% mals_z=zscore(mals);
% 
% owBoldTs = regressOut(iwBoldTs.',mals_z.',1).';
% myOutputBold=ts2vol(owBoldTs,logical(ones(nx,ny,nz)));
% myoBoldTs=vol2ts(myOutputBold,roiMask);
% 
% %% ADD FILTER AT 0.1 Hz
% fs=1/TR;
% [bf,af]=butter(4,0.1/(fs/2));
% 
% nPad=30;
% tmp=cat(1,repmat(myoBoldTs(1,:),nPad,1),myoBoldTs);
% tmpOut=filter(bf,af,tmp,[],1);
% myfoBoldTs=tmpOut(nPad+1:end,:);
% 
% figure;
% subplot(211);
% plot(mean(oBoldTs,2));
% subplot(212);
% plot(mean(myfoBoldTs,2));
% corrcoef(mean(oBoldTs,2),mean(myfoBoldTs,2))

% 
% %% 
% % segment skull stripped anatomical
% str=['!3dSeg -anat ' basePath 'anat_ns+orig -mask AUTO -classes ''CSF ; GM ; WM'' -bias_classes ''GM ; WM'' -bias_fwhm 25 -mixfrac UNI -main_N 5 -blur_meth BFT -prefix ' basePath];
% eval(str);
% 
% %%
% % coregister segmentation with talairach
% origPath=pwd;
% cd(basePath);
% str='!@auto_tlrc -apar anat+tlrc -input Classes+orig.BRIK';
% eval(str);
% cd(origPath);
% 
% %%
% % resample talairach segmentation on bold grid
% origPath=pwd;
% cd(basePath);
% str='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_Classes -input Classes+tlrc';
% eval(str);
% cd(origPath);
% 
% 
% 
% 
% 
% %%
% [~, outputBold, Info, ~] = BrikLoad (outputBrik);
% 
% %%

% 
% % % draw

% 
% 
% %%
% % regress out CSF
% [~, segBold, Info, ~] = BrikLoad (segBrik);
% 
% CSF=1;
% csfMask=segBold==CSF;  
% csfBoldTs=vol2ts(outputBold,csfMask);
% [Uc,Sc,Vc]=svd(csfBoldTs,0);
% 
% %%
% WM=3;  
% wmMask=segBold==WM;
% wmBoldTs=vol2ts(outputBold,wmMask);
% [Uw,Sw,Vw]=svd(wmBoldTs,0);
% 
% %%
% wbBoldTs=vol2ts(outputBold);
% ooBoldTs=regressOut(oBoldTs.',Uc(:,1:5).',1).';
% oooBoldTs=regressOut(oBoldTs.',Uw(:,1:5).',1).';
% %%
% % 
% 
% %% draw
% figure;
% subplot(211);
% plot(mean(oBoldTs,2));
% subplot(212);
% %plot(mean(ooBoldTs,2));
% plot(mean(oooBoldTs,2));
% 
% %%
% figure; hold on
% plot(mean(csfBoldTs,2));
% plot(mean(wmBoldTs,2));
% 
% 
% %precomputedFilename=['../data/precomputed/allBoldsRoi 22-Apr-2018'];
% %load(precomputedFilename,'allBoldsRoiOut','allOnsets','TR','TEs','subjStr','nSubjects','nEchos');
% % sBold=mean(squeeze(allBoldsRoiOut{sIndx}(1,:,:)),2);
% % % pull in motion alignment parameters
% % malFile='../data/S02/NII/dsbold_e1_vr_motion.1D';
% % fid=fopen(malFile);
% % mals=textscan(fid,'%f %f %f %f %f %f');
% % mals=cell2mat(mals);
% % dM=cat(1,zeros(1,6),diff(mals));
% % figure
% % subplot(211);
% % plot(sBold);
% % subplot(212);
% % plot(sum(dM.^2,2));