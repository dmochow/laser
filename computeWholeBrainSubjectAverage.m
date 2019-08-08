% May 4th, 2018
% setting up for the whole-brain analysis
%
% take the average of everyone's bold 
% this script after preprocess bold
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

subjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19'};
nSubjects=numel(subjStr);

brainMaskFilename='resampled_brain_mask+tlrc';
epiFilename_1='smnudsbold_e1_tlrc_al+tlrc';
epiFilename_2='smnudsbold_e2_tlrc_al+tlrc';
epiFilename_3='smnudsbold_e3_tlrc_al+tlrc';
biopacFilename='biopac.mat';
nEchos=3;
TR=2.8;
nTR=645;
nx=64;ny=76;nz=60; % BOLD dimensions
nominalOnsetTR=214; % 

muBold_1=zeros(nx,ny,nz,nTR);
muBold_2=zeros(nx,ny,nz,nTR);
muBold_3=zeros(nx,ny,nz,nTR);

pathToWriteData='../data/SAVG/NII/';
%%
for s=1:nSubjects
    s
    path=['../data/' subjStr{s} '/NII/'];
    biopacPath=['../data/' subjStr{s} '/BIOPAC/'];
      
    %%
    % echo 1
    [~, brainMask, Info, ~] = BrikLoad (fullfile(path,brainMaskFilename));
    brainMask=logical(brainMask);

    % load in bold and make 2D series
    [~, bold_1, Info, ~] = BrikLoad (fullfile(path,epiFilename_1));
    bold1_2D=vol2ts(bold_1,brainMask);
    
    % global signal regression  
    muBold1_2D=mean(bold1_2D,2);
    bold1_2D_r=regressOut(bold1_2D.',muBold1_2D.',1).';  
    bold1_vol=ts2vol(bold1_2D_r,brainMask);
    
    %%
    % echo 2
    [~, bold_2, Info, ~] = BrikLoad (fullfile(path,epiFilename_2));
    bold2_2D=vol2ts(bold_2,brainMask);
    muBold2_2D=mean(bold2_2D,2);
    bold2_2D_r=regressOut(bold2_2D.',muBold2_2D.',1).';  
    bold2_vol=ts2vol(bold2_2D_r,brainMask);
    
    %%
    % echo 3
    [~, bold_3, Info, ~] = BrikLoad (fullfile(path,epiFilename_3));
    bold3_2D=vol2ts(bold_3,brainMask);
    muBold3_2D=mean(bold3_2D,2);
    bold3_2D_r=regressOut(bold3_2D.',muBold3_2D.',1).';  
    bold3_vol=ts2vol(bold3_2D_r,brainMask);
    
    % laser onset times from biopac
    load(fullfile(biopacPath,biopacFilename),'data');
    
    [~,maxindsample]=max(diff(data(:,1)));
    onsetTimeSec=maxindsample/1000;
    onsetTimeTR=round(onsetTimeSec/TR);
        
    shiftVal=nominalOnsetTR-onsetTimeTR;
    muBold_1=muBold_1+circshift(bold1_vol,[1 1 1 shiftVal])/nSubjects;
    muBold_2=muBold_2+circshift(bold2_vol,[1 1 1 shiftVal])/nSubjects;
    muBold_3=muBold_3+circshift(bold3_vol,[1 1 1 shiftVal])/nSubjects;
    
end

%%
[err1,ErrMessage1,~]=WriteBrikWrap(pathToWriteData,muBold_1,Info,'subjectMeanBold_e1','tlrc');
[err2,ErrMessage2,~]=WriteBrikWrap(pathToWriteData,muBold_2,Info,'subjectMeanBold_e2','tlrc');
[err3,ErrMessage3,~]=WriteBrikWrap(pathToWriteData,muBold_3,Info,'subjectMeanBold_e3','tlrc');

%precomputedFilename=['../data/precomputed/allBoldsRoi ' date];
%save(precomputedFilename,'allBoldsRoi','allBoldsGlobal','allOnsets','allBoldsControlRoi','TR','TEs','subjStr','nSubjects','nEchos');

