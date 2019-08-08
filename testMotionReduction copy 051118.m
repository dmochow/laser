clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%%
subjIndx='S11';
basePath=['../data/' subjIndx '/NII/'];
inputBrik=[basePath 'dsbold_e1_tlrc_al+tlrc.BRIK'];
outputBrik=[basePath 'mynudsbold_e1_tlrc_al+tlrc.BRIK'];
segBrik=[basePath 'resampled_Classes+tlrc.BRIK'];
roiMaskFilename=[basePath 'muroi+tlrc.BRIK'];
anatFilename=[basePath 'anat_ns+orig.BRIK'];
motionFile=[basePath 'dsbold_e1_vr_motion.1D'];
[~, inputBold, Info, ~] = BrikLoad (inputBrik);
[~, roiMask, Info, ~] = BrikLoad (roiMaskFilename);
roiMask=logical(roiMask);

%% 
% segment skull stripped anatomical
str=['!3dSeg -anat ' basePath 'anat_ns+orig -mask AUTO -classes ''CSF ; GM ; WM'' -bias_classes ''GM ; WM'' -bias_fwhm 25 -mixfrac UNI -main_N 5 -blur_meth BFT -prefix ' basePath];
eval(str);

%%
% coregister segmentation with talairach
origPath=pwd;
cd(basePath);
str='!@auto_tlrc -apar anat+tlrc -input Classes+orig.BRIK';
eval(str);
cd(origPath);

%%
% resample talairach segmentation on bold grid
origPath=pwd;
cd(basePath);
str='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_Classes -input Classes+tlrc';
eval(str);
cd(origPath);


%%
% regress out realignment parameters
str=['delete ' basePath 'mynudsbold_e1_tlrc_al+tlrc.BRIK'];
eval(str);
str=['delete ' basePath 'mynudsbold_e1_tlrc_al+tlrc.HEAD'];
eval(str);

str=['!3dTproject -input ' inputBrik ' -prefix ' basePath 'mynudsbold_e1_tlrc_al -ort ' motionFile];
eval(str);


%%
[~, outputBold, Info, ~] = BrikLoad (outputBrik);

%%
% compare input to output
iBoldTs=vol2ts(inputBold,roiMask);
oBoldTs=vol2ts(outputBold,roiMask);

% % draw
% figure;
% subplot(211);
% plot(mean(iBoldTs,2));
% subplot(212);
% plot(mean(oBoldTs,2));


%%
% regress out CSF
[~, segBold, Info, ~] = BrikLoad (segBrik);

CSF=1;
csfMask=segBold==CSF;  
csfBoldTs=vol2ts(outputBold,csfMask);
[Uc,Sc,Vc]=svd(csfBoldTs,0);

%%
WM=3;  
wmMask=segBold==WM;
wmBoldTs=vol2ts(outputBold,wmMask);
[Uw,Sw,Vw]=svd(wmBoldTs,0);

%%
wbBoldTs=vol2ts(outputBold);
ooBoldTs=regressOut(oBoldTs.',Uc(:,1:5).',1).';
oooBoldTs=regressOut(oBoldTs.',Uw(:,1:5).',1).';
%%
% 

%% draw
figure;
subplot(211);
plot(mean(oBoldTs,2));
subplot(212);
%plot(mean(ooBoldTs,2));
plot(mean(oooBoldTs,2));

%%
figure; hold on
plot(mean(csfBoldTs,2));
plot(mean(wmBoldTs,2));


%precomputedFilename=['../data/precomputed/allBoldsRoi 22-Apr-2018'];
%load(precomputedFilename,'allBoldsRoiOut','allOnsets','TR','TEs','subjStr','nSubjects','nEchos');
% sBold=mean(squeeze(allBoldsRoiOut{sIndx}(1,:,:)),2);
% % pull in motion alignment parameters
% malFile='../data/S02/NII/dsbold_e1_vr_motion.1D';
% fid=fopen(malFile);
% mals=textscan(fid,'%f %f %f %f %f %f');
% mals=cell2mat(mals);
% dM=cat(1,zeros(1,6),diff(mals));
% figure
% subplot(211);
% plot(sBold);
% subplot(212);
% plot(sum(dM.^2,2));