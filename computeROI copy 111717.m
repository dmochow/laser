clear all; close all; clc
addpath(genpath('../../NIfTI_20140122'));

%
niiFilename='/Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat.nii';

%
nii=load_untouch_nii(niiFilename);
img=nii.img;
% dim 1: sagittal
% dim 2: coronal
% dim 3: axial
% figure;
% imagesc(squeeze(img(:,128,:)));
% colormap bone

%%
nx=size(img,1);
ny=size(img,2);
nz=size(img,3);
[X,Y,Z]=ndgrid(1:nx,1:ny,1:nz);

%% left
llf=[130 10 163]+[1 1 1]; % adding 1s because AFNI indexes from 0
ulf=[130 11 171]+[1 1 1];
ulb=[133 1 174]+[1 1 1];
isRight = makeRoiMask(ulf,llf,ulb,X,Y,Z,1);

%% right
lrf=[161 20 161]+[1 1 1]; % adding 1s because AFNI indexes from 0
urf=[161 22 168]+[1 1 1];
urb=[167 1 172]+[1 1 1];
isLeft = makeRoiMask(lrf,urf,urb,X,Y,Z,-1);

%% define "below" mask
x1=ulf;
x2=ulb;
x3=urb;
y1=x1; y1(2)=256-y1(2);
y2=x2; y2(2)=256-y2(2);
y3=x3; y3(2)=256-y3(2);
cross_u=cross(y1-y2,y1-y3);
cross_u=cross_u/norm(cross_u);
if cross_u(3)>0 % pointing up
    deld=20;
else
    deld=-20;     % pointing down
end
z1=y1+deld*cross_u;
z2=y2+deld*cross_u;
z3=y3+deld*cross_u;
q1=z1; q1(2)=256-z1(2);
q2=z2; q2(2)=256-z2(2);
q3=z3; q3(2)=256-z3(2);
isBelow = makeRoiMask(q1,q2,q3,X,Y,Z,-1);

% figure; 
% subplot(121);hold on
% plot3(y1(1),y1(2),y1(3),'*b');
% plot3(y2(1),y2(2),y2(3),'*b');
% plot3(y3(1),y3(2),y3(3),'*b');
% plot3(z1(1),z1(2),z1(3),'*r');
% plot3(z2(1),z2(2),z2(3),'*r');
% plot3(z3(1),z3(2),z3(3),'*r');
% xlim([1 224]); ylim([1 256]); zlim([1 256]);
% view([45 18]);
% 
% subplot(122); hold on
% plot3(x1(1),x1(2),x1(3),'*b');
% plot3(x2(1),x2(2),x2(3),'*b');
% plot3(x3(1),x3(2),x3(3),'*b');
% plot3(q1(1),q1(2),q1(3),'*r');
% plot3(q2(1),q2(2),q2(3),'*r');
% plot3(q3(1),q3(2),q3(3),'*r');
% xlim([1 224]); ylim([1 256]); zlim([1 256]);
% view([45 18]);


%% define "above" mask
% here derived as parallel to the below mask
x1=ulf;
x2=ulb;
x3=urb;
y1=x1; y1(2)=256-y1(2);
y2=x2; y2(2)=256-y2(2);
y3=x3; y3(2)=256-y3(2);
cross_u=cross(y1-y2,y1-y3);
cross_u=cross_u/norm(cross_u);
if cross_u(3)>0 % pointing up
    deld=-20;
else
    deld=20;     % pointing down
end
z1=y1+deld*cross_u;
z2=y2+deld*cross_u;
z3=y3+deld*cross_u;
q1=z1; q1(2)=256-z1(2);
q2=z2; q2(2)=256-z2(2);
q3=z3; q3(2)=256-z3(2);
isAbove = makeRoiMask(q1,q2,q3,X,Y,Z,0);


%% make mask
mask=zeros(size(img));
mask(and(and(and(isLeft,isRight),isAbove),isBelow))=255;

% write mask to nifti
nii2=nii;
nii2.img=mask;
nii2filename='/Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/roiMask.nii';
save_untouch_nii(nii2,nii2filename);

%%
figure;
subplot(211);
imagesc(img(:,:,174));
colormap bone
subplot(212);
%imagesc(isLeftImg(:,:,170));
imagesc(mask(:,:,256));
colormap bone


