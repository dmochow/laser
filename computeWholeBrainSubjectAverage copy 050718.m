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
      
    [~, brainMask, Info, ~] = BrikLoad (fullfile(path,brainMaskFilename));
    [~, bold_1, Info, ~] = BrikLoad (fullfile(path,epiFilename_1));
    [~, bold_2, Info, ~] = BrikLoad (fullfile(path,epiFilename_2));
    [~, bold_3, Info, ~] = BrikLoad (fullfile(path,epiFilename_3));
    
    %bold1_2D=permute(reshape(bold_1,[nx*ny*nz nTR]),[2 1]);
    bold1_2D=permute(bold_1,[4 1 2 3]); bold1_2D=bold1_2D(:,:);
    bold1_2D_masked=bold1_2D(:,find(brainMask));
    
    bold2_2D=permute(reshape(bold_2,[nx*ny*nz nTR]),[2 1]);
    bold2_2D_masked=bold2_2D(:,find(brainMask));
    
    bold3_2D=permute(reshape(bold_3,[nx*ny*nz nTR]),[2 1]);
    bold3_2D_masked=bold3_2D(:,find(brainMask));
    
    % global signal regression  
    muBold1_2D_masked=mean(bold1_2D_masked,2);
    bold1_2D_masked_r=regressOut(bold1_2D_masked.',muBold1_2D_masked.',1).';  
    
    muBold2_2D_masked=mean(bold2_2D_masked,2);
    bold2_2D_masked_r=regressOut(bold2_2D_masked.',muBold2_2D_masked.',1).';  
    
    muBold3_2D_masked=mean(bold3_2D_masked,2);
    bold3_2D_masked_r=regressOut(bold3_2D_masked.',muBold3_2D_masked.',1).'; 

    %%
    % unmask and put back to 4D
    tBold1=bold_1;
    tBold1_2D=permute(tBold1,[4 1 2 3]); tBold1_2D=tBold1_2D(:,:);
    tBold1_2D(:,find(brainMask))=bold1_2D_masked;
    tmp=reshape(tBold1,[nTR nx ny nz]);
    tmp2=permute(tmp,[2 3 4 1]);
    sum(abs(tmp2(:)-bold_1(:)))
    
    %%
    tmp=permute(bold1_2D_masked,[2 1]);
    for n=1:nTR
        tBold1=zeros(nx*ny*nz,nTR);
        tBold1(find(brainMask),n)=tmp(:,n);
    end
    tBold1=reshape(tBold1,[nx ny nz nTR]);
    %%

    
    % laser onset times from biopac
    load(fullfile(biopacPath,biopacFilename),'data');
    
    [~,maxindsample]=max(diff(data(:,1)));
    onsetTimeSec=maxindsample/1000;
    onsetTimeTR=round(onsetTimeSec/TR);
        
    shiftVal=nominalOnsetTR-onsetTimeTR;
    muBold_1=muBold_1+circshift(bold_1,[1 1 1 shiftVal])/nSubjects;
    muBold_2=muBold_2+circshift(bold_2,[1 1 1 shiftVal])/nSubjects;
    muBold_3=muBold_3+circshift(bold_3,[1 1 1 shiftVal])/nSubjects;
    
end

%%
[err1,ErrMessage1,~]=WriteBrikWrap(pathToWriteData,muBold_1,Info,'subjectMeanBold_e1','tlrc');
[err2,ErrMessage2,~]=WriteBrikWrap(pathToWriteData,muBold_2,Info,'subjectMeanBold_e2','tlrc');
[err3,ErrMessage3,~]=WriteBrikWrap(pathToWriteData,muBold_3,Info,'subjectMeanBold_e3','tlrc');

%precomputedFilename=['../data/precomputed/allBoldsRoi ' date];
%save(precomputedFilename,'allBoldsRoi','allBoldsGlobal','allOnsets','allBoldsControlRoi','TR','TEs','subjStr','nSubjects','nEchos');

