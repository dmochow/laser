% analyze whole brain average

clear all; close all; clc

path='../data/SAVG/NII/';
epiFilename_1='subjectMeanBold_e1+tlrc.BRIK';
epiFilename_2='subjectMeanBold_e2+tlrc.BRIK';
epiFilename_3='subjectMeanBold_e3+tlrc.BRIK';
brainMaskFilename='muRoi+tlrc.BRIK';

%%
[~, brainMask, Info, ~] = BrikLoad (fullfile(path,brainMaskFilename));
brainMask=logical(brainMask);

%%
% echo 1
[~, bold_1, Info, ~] = BrikLoad (fullfile(path,epiFilename_1));
bold1_2D=vol2ts(bold_1,brainMask);
bold1_2D_=vol2ts(bold_1); muBold1_2D=mean(bold1_2D_,2);
bold1_2D=regressOut(bold1_2D.',muBold1_2D.',1).';

%%
% echo 2
[~, bold_2, Info, ~] = BrikLoad (fullfile(path,epiFilename_2));
bold2_2D=vol2ts(bold_2,brainMask);
bold2_2D_=vol2ts(bold_2); muBold2_2D=mean(bold2_2D_,2);
bold2_2D=regressOut(bold2_2D.',muBold2_2D.',1).';

%%
% echo 3
[~, bold_3, Info, ~] = BrikLoad (fullfile(path,epiFilename_3));
bold3_2D=vol2ts(bold_3,brainMask);
bold3_2D_=vol2ts(bold_3); muBold3_2D=mean(bold3_2D_,2);
bold3_2D=regressOut(bold3_2D.',muBold3_2D.',1).';

%%
% echo mean
muBold_2D=(bold1_2D+bold2_2D+bold3_2D)/3;

%bold1_2D_c=vol2ts(bold_1,~brainMask);
figure; 
subplot(221); plot(mean(bold1_2D,2));
subplot(222); plot(mean(bold2_2D,2));
subplot(223); plot(mean(bold3_2D,2));
subplot(224); plot(mean(muBold_2D,2));
