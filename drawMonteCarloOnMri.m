clear all; close all; clc

%%
% data to project
load ../data/mcml/Arz.mat
x=[r(1) r(end)];
y=[-z(end) z(end)];
Adraw=cat(2,Abs(end:-1:1,:).',Abs.');

%%
laserOrigin=load('../data/S11/NII/laserOrigin.mat');

%%
anatFilename='../data/S11/NII/anat+orig';
[err, V, InfoMask, ErrMessage] = BrikLoad(anatFilename);

sl=[100 100 100];

figure;
subplot(221)
imagesc(squeeze(V(sl(1),:,:))); colormap bone
subplot(222)
imagesc(squeeze(V(:,sl(2),:))); colormap bone
subplot(223)
imagesc(squeeze(V(:,:,sl(3)))); colormap bone
subplot(224);
imagesc(y,x,Adraw); colormap jet