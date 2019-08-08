clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

% point to the data
path='../data/S10/NII/';
maskFilename='resampled_roi+tlrc';
anatFilename='resampled_anat+tlrc';
epiFilename_1='bold_e1_tlrc_al+tlrc';
epiFilename_2='bold_e2_tlrc_al+tlrc';
epiFilename_3='bold_e3_tlrc_al+tlrc';
TEs=[13 34 55];

% read in the data
[err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,maskFilename));
[err, anat, Info, ErrMessage] = BrikLoad (fullfile(path,anatFilename));
[err, epi_1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
[err, epi_2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
[err, epi_3, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_3));

% compute the mask
finalMask=mask>0 & anat>0; % in the laser beam AND in the brain

% grab the bold in the ROI
tmp=permute(epi_1,[4 1 2 3]);
epi(1,:,:)=tmp(:,:);

tmp=permute(epi_2,[4 1 2 3]);
epi(2,:,:)=tmp(:,:);

tmp=permute(epi_3,[4 1 2 3]);
epi(3,:,:)=tmp(:,:);

%%
epi_roi=epi(:,:,finalMask);

%%

% perform SVD on single echo
epi1=squeeze(epi_roi(1,:,:)).';
epi2=squeeze(epi_roi(2,:,:)).';
epi3=squeeze(epi_roi(3,:,:)).';
[U,S,V]=svd(epi2,0);
% U are the spatial filters

% now project the spatial filters onto each echo
B1=U'*epi1; % for each component, the expression at each time point
B2=U'*epi2;
B3=U'*epi3;
% NB: component X time 

%% compute F-scores

%% segregate into BOLD and non-BOLD


%%
nComps2Show=5;
figure
for c=1:nComps2Show
    subplot(nComps2Show,2,(c-1)*2+1);
    plot(U(:,c));
    subplot(nComps2Show,2,c*2);
    plot(V(:,c));
end    
    

%% 
% plot TE dependence of components
cIndx=6;
figure;
plot(TEs,[B1(cIndx,:); B2(cIndx,:); B3(cIndx,:)]);

%%







% look at time courses
%%
figure;
plot(V(:,5));


% look at spectrum


% %%
% nReg=2;
% nComp=2;
% [dataOut,W,A,Rxx,Ryy,Rxy,dGen] = rcaRun(data,nReg,nComp,[],[],0,[]);
% 
% 
% %epi2Dmasked=epi2D(:,mask1D>0);
% % 
% % %%
% % meanRoiTs=mean(epi2Dmasked,2);
% % figure;
% % plot(meanRoiTs);
% % 
% % %%
% % [U,S,V]=svd(epi2Dmasked);
% % figure;
% % for c=1:12
% %     subplot(3,4,c);
% %     plot(U(:,c),'k');
% % end

%% assume biopac in memory as data
% [~,maxindsample]=max(diff(data(:,1)));
% onsetTimeSec=maxindsample/1000;
% onsetTimeTR=onsetTimeSec/2.8;

