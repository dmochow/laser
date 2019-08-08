clear all; close all; clc
% figure out whether 3dtproject scales the data or if it's inherent to the
% regression
subjStr='S13'; % our test subject
path=['../data/' subjStr '/NII/'];
epiFilename_1='dsbold_e3_tlrc_al+tlrc';
epiFilename_2='nudsbold_e3_tlrc_al+tlrc';
nuFilename='dsbold_e3_vr_motion.1D';
[err, data1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
[err, data2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
%%
% grab nuisance regressors
Y = importdata(fullfile(path,nuFilename));
Y2=Y-repmat(mean(Y,1),size(Y,1),1);
Y3=Y2./repmat(std(Y2,[],1),size(Y2,1),1);

[U,S,V]=svd(Y3,0);
%%
% project out
voxelIn=squeeze(data1(36-1,75-10-1,41-1,:));
voxelOut=squeeze(data2(36-1,75-10-1,41-1,:));
%myvoxelOut = (eye(size(Y3,1))-Y3*pinv(Y3))*voxelIn;
Q=eye(size(U,1))-Y*pinv(Y);
%myvoxelOut = Q*voxelIn;
myvoxelOut=regressOut(voxelIn,Y3);
%myvoxelOut=voxelIn-Y3*pinv(Y3)*voxelIn;

figure;
subplot(311)
plot(voxelIn);
subplot(312)
plot(voxelOut);
subplot(313)
plot(myvoxelOut);

% %%
% nfft=1024;
% fs=1/2.8;
% freqs=(0:nfft-1)/nfft*fs;
% Vin=fft(voxelIn,nfft);
% Vout=fft(voxelOut,nfft);
% myVout=fft(myvoxelOut,nfft);
% figure;
% subplot(311)
% plot(freqs,abs(Vin));  xlim([0 fs/2]);
% subplot(312);
% plot(freqs,abs(Vout)); xlim([0 fs/2]);
% subplot(313);
% plot(freqs,abs(myVout)); xlim([0 fs/2]);

