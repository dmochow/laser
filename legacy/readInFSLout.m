clear all; close all; clc
% read in filtered data from FSL and try to analyze it
% high-pass filtering
% experimental design
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));
addpath(genpath('./vol3d'));

%%
subjIndx=2;
alpha=0.05; % for output nii
nPerms=500; % for significance testing
%%
switch subjIndx
    case 1
        maskFilename='../data/S01/FSL.ica/filtered_func_data.ica/mask.nii';
        dataFilename='../data/S01/FSL.ica/filtered_func_data.nii';
        tOn=401:800;
        tOff=[1:400 801:1200];
        fs=1/1.5;  % 
        outFilename='../data/S01/FSL.ica/signalRatio_whighpass.nii';
        outFilename2='../data/S01/FSL.ica/pvalues_whighpass.nii';
    case 2
        maskFilename='../data/S02/FSL.ica/filtered_func_data.ica/mask.nii';
        dataFilename='../data/S02/FSL.ica/filtered_func_data.nii';
        tOn=216:430;
        tOff=[1:215 431:645];
        fs=1/2.8;  
        outFilename='../data/S02/FSL.ica/signalRatio_whighpass.nii';
        outFilename2='../data/S02/FSL.ica/pvalues_whighpass.nii';
end

%%
maskNii=load_untouch_nii(maskFilename);
mask=logical(maskNii.img);

%%
dataNii=load_untouch_nii(dataFilename);
data=double(dataNii.img);
data=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
brainData=data(mask(:),:);

%%
% add high-pass filtering
% [filt_b,filt_a]=butter(6,0.01/fs*2,'high');
% brainDataFilt=filter(filt_b,filt_a,brainData,[],2);
%% 
% for each voxel in the brain, compute a ratio of the signal during stim to
% the signal before/after stim
% signalOn=mean(brainData(:,tOn),2);
% signalOff=mean(brainData(:,tOff),2);
signalOn=std(brainData(:,tOn),[],2);
signalOff=std(brainData(:,tOff),[],2);
signalRatio=signalOn./signalOff;

%%
% basic figure
clrs=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250];
figure
[counts,bins]=hist(signalRatio,1000);
bar(bins,counts/sum(counts));
xlabel('Laser On : Laser Off');
ylabel('Relative Frequency');


%%
% permutation testing
mock_signalRatio=zeros(size(signalRatio,1),nPerms);
nOn=numel(tOn);
nOff=numel(tOff);
for p=1:nPerms
    p
%     rand_design = zeros(nOn+nOff,1);
%     rand_design(randperm(nOn+nOff, nOn)) = 1;
%     mock_tOn=find(rand_design);
%     mock_tOff=find(~rand_design);
% 
%     mock_signalOn=std(brainData(:,mock_tOn),[],2);
%     mock_signalOff=std(brainData(:,mock_tOff),[],2);
%     mock_signalRatio(:,p)=mock_signalOn./mock_signalOff;
    
    mock_brainData=surrogateResponseGenerator(brainData')';
    
    mock_signalOn=std(mock_brainData(:,tOn),[],2);
    mock_signalOff=std(brainData(:,tOff),[],2);
    mock_signalRatio(:,p)=mock_signalOn./mock_signalOff;
    
end

%% 
% compute pvals
nVoxels=size(signalRatio,1);
pvals=zeros(nVoxels,1);
for v=1:nVoxels
    pvals(v)=mean(mock_signalRatio(v,:)>signalRatio(v));
end

%%
% correct p-values for multiple comparisons
[p_fdr, p_masked] = fdr( pvals, alpha);


%%
tmp=zeros(size(mask(:)));
tmp(mask(:))=signalRatio;
ratioVolume=reshape(tmp,[size(mask,1),size(mask,2),size(mask,3)]);

tmp=zeros(size(mask(:)));
tmp(mask(:))=p_masked;
pvalVolume=reshape(tmp,[size(mask,1),size(mask,2),size(mask,3)]);



% ratioVolume_withNan=ratioVolume;
% ratioVolume_withNan(~mask)=NaN;
% ratiosLeft=ratioVolume_withNan(1:45,:,:); nanmean(ratiosLeft(:))
% ratiosRight=ratioVolume_withNan(46:end,:,:);nanmean(ratiosRight(:))
% ratiosFront=ratioVolume_withNan(:,1:45,:); nanmean(ratiosFront(:))
% ratiosBack=ratioVolume_withNan(:,46:end,:);nanmean(ratiosBack(:))
% ratiosUp=ratioVolume_withNan(:,:,31:60); nanmean(ratiosUp(:))
% ratiosDown=ratioVolume_withNan(:,:,1:30);nanmean(ratiosDown(:))

%%
figure;
subplot(121); vol3d('cdata',ratioVolume);
subplot(122); vol3d('cdata',pvalVolume);
% axis tight;  daspect([1 1 .4])
% alphamap('rampup');
% alphamap(.5 .* alphamap);

%%
% write output nii
outNii=maskNii;
outNii.img=ratioVolume;
save_untouch_nii(outNii,outFilename);

outNii2=maskNii;
outNii2.img=pvalVolume;
save_untouch_nii(outNii2,outFilename2);


%%
% figure;
% plot(mean(brainData.^2,1))