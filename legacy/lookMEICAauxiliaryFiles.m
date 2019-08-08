clear all; close all; clc

addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));
subjIndx=4;

%/Users/jacek/Documents/MATLAB/LLLT/data/S04/NII/mebold2go/meica.bold_e123/TED

maskFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/meica.bold_e123/eBvrmask.nii'];
dataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn.nii'];
dntsFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/meica.bold_e123/TED/dn_ts_OC.nii'];
biopacFilename=['../data/S0' num2str(subjIndx) '/BIOPAC/greg.mat'];
tmp=['../data/S0' num2str(subjIndx) '/NII/mebold2go/meica.bold_e123/TED/hik_ts_OC.nii'];
outFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn_pvals.nii'];
saveDataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/matData'];
TR=2.8;
fs=1/2.8;
DUR_frames=645; % duration of BOLD recording in frames

%%
maskNii=load_untouch_nii(maskFilename);
mask=logical(maskNii.img);

%%
dataNii=load_untouch_nii(dntsFilename);
data=double(dataNii.img);
data_4D=data;
data=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
brainData=data(mask(:),:);
