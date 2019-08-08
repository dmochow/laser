clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

% this file makes a mock ROI that is contralateral to the actual laser!


subjStr='S11';

% important parameters
SCALP_INT_THRESH=150; % for finding edge of scalp
MAXITER=1000;
STEPSIZE=0.1; % 100 microns
RADMASK=19.2; % lateral laser reach
DEPTHMASK=26.5; % depth laser reach

switch subjStr
    case 'S04'
        urf=[158 255-26 164]+[1 1 1];
        ulf=[123 255-24 166]+[1 1 1];
        urb=[164 255-8 164]+[1 1 1];
        lrf=[156 255-28 151]+[1 1 1];
        
    case 'S05'
        urf=[150 255-22 159]+[1 1 1];
        ulf=[116 255-20 161]+[1 1 1];
        urb=[156 255-1 164]+[1 1 1];
        lrf=[150 255-22 151]+[1 1 1];
        
    case 'S06' 
        urf=[161 255-20 170]+[1 1 1];
        ulf=[128 255-17 170]+[1 1 1];
        urb=[169 255-2 170]+[1 1 1];
        lrf=[161 255-19 162]+[1 1 1];
        
    case 'S07'
        urf=[158 255-27 165]+[1 1 1];
        ulf=[126 255-20 171]+[1 1 1];
        urb=[165 255-6 167]+[1 1 1];
        lrf=[156 255-27 156]+[1 1 1];
        
    case 'S09'
        ulf=[126 255-43 179]+[1 1 1];
        urf=[147 255-49 179]+[1 1 1];
        urb=[154 255-35 192]+[1 1 1];
        lrf=[138 255-17 149]+[1 1 1];
        
    case 'S10'
        ulf=[141 255-23 181]+[1 1 1];
        urf=[158 255-28 177]+[1 1 1];
        urb=[165 255-9 180]+[1 1 1];
        lrf=[155 255-18 137]+[1 1 1];
        
    case 'S11'
        midsag=111;
        ulf=[2*midsag-127 255-48 197]+[1 1 1];
        urf=[2*midsag-148 255-55 193]+[1 1 1];
        urb=[2*midsag-161 255-25 213]+[1 1 1];
        lrf=[2*midsag-149 255-30 158]+[1 1 1];
        
    case 'S12'
        urf=[151 255-47 198]+[1 1 1]; %
        ulf=[129 255-41 202]+[1 1 1]; %
        lrf=[151 255-34 159]+[1 1 1]; %
        urb=[152 255-32 203]+[1 1 1]; %
        
    case 'S13'
        urf=[146 255-31 182]+[1 1 1]; %
        ulf=[123 255-28 186]+[1 1 1]; %
        urb=[153 255-3 192]+[1 1 1];
        lrf=[144 255-19 143]+[1 1 1]; %      
end
pathToData=['../data/' subjStr '/NII/'];

% % try to compute brain mask in afni
% origPath=pwd;
% cd(pathToData);
% if exist('brain_mask+orig.BRIK','file')
%     delete brain_mask+orig.BRIK;
%     delete brain_mask+orig.HEAD;
% end
% str1='!3dSkullStrip anat.nii';
% str2='!3dAutomask -prefix brain_mask skull_strip_out+orig';
% eval(str1);
% eval(str2);
% cd(origPath);

% processing begins here
pathToData=['../data/' subjStr '/NII/'];
niiFilename=['../data/' subjStr '/NII/anat.nii'];
anatFilename=['../data/' subjStr '/NII/anat+orig'];
% brainMaskFilename=['../data/' subjStr '/NII/brain_mask+orig'];
outNiifilename=['../data/' subjStr '/NII/controlRoiMask.nii'];

laserOrigin = getLaserOrigin(ulf,urf,urb); % origin in afni space, where second dimension is A->P

%%
%project towards the skin until we hit it
[err, anat, Info, ErrMessage] = BrikLoad (anatFilename);

% define front face
crossprod=cross((ulf-urf),(ulf-lrf));
crossprod=crossprod/norm(crossprod);
if crossprod(2)>0 % if pointing anterior
    crossprod=-crossprod;
end
xo=laserOrigin;
imageInt=zeros(MAXITER,1);
for i=1:MAXITER
    xo=xo+crossprod*STEPSIZE;
    rxo=round(xo);
    try
        imageInt(i)=anat(rxo(1),rxo(2),rxo(3));
    catch
        imageInt(i)=imageInt(i-1);
    end
end

dImageInt=[0;diff(imageInt)];

% figure;
% subplot(211);
% plot(1:nIter,imageInt);
% x=ginput(1);
% maxind=x(1);
% subplot(212);
% plot(1:nIter,dImageInt);
%
% %laserOrigin=laserOrigin+maxind*ss*crossprod

findind=find(dImageInt>SCALP_INT_THRESH,1);
laserOrigin=laserOrigin+findind*STEPSIZE*crossprod

save(['../data/' subjStr '/NII/controlLaserOrigin.mat'],'laserOrigin','findind','imageInt');

%%
%[err, mask, Info, ErrMessage] = BrikLoad (brainMaskFilename);
%mridim=size(mask);
mridim=size(anat);
[x,y,z]=ndgrid(0:mridim(1)-1,0:mridim(2)-1,0:mridim(3)-1);
LOx=laserOrigin(1)*ones(size(x));
LOy=laserOrigin(2)*ones(size(y));
LOz=laserOrigin(3)*ones(size(z));

%%
% translation
T=eye(4);
T(1,4)=-laserOrigin(1);
T(2,4)=-laserOrigin(2);
T(3,4)=-laserOrigin(3);

xt=T(1,1)*x+T(1,2)*y+T(1,3)*z+T(1,4)*1;
yt=T(2,1)*x+T(2,2)*y+T(2,3)*z+T(2,4)*1;
zt=T(3,1)*x+T(3,2)*y+T(3,3)*z+T(3,4)*1;

%%
% rotation
u=(urf-ulf)/norm(urf-ulf);  % LR
v=(urf-urb)/norm(urf-urb);  % PA
w=(urf-lrf)/norm(urf-lrf);  % IS
R=eye(4);
R(1:3,1)=u;
R(1:3,2)=v;
R(1:3,3)=w;
R=inv(R);

xr=R(1,1)*xt+R(1,2)*yt+R(1,3)*zt+R(1,4)*1;
yr=R(2,1)*xt+R(2,2)*yt+R(2,3)*zt+R(2,4)*1;
zr=R(3,1)*xt+R(3,2)*yt+R(3,3)*zt+R(3,4)*1;


%%
roiMask=zeros(size(anat));
roiMask( (xr.^2+zr.^2)<RADMASK.^2 & yr>0 & yr<DEPTHMASK )=1;

%% write out ROI mask
InfoWrite=Info;
opt.Scale=0;
opt.Prefix='controlRoi';
opt.View='orig';
opt.Verbose=1;
opt.AppendHistory=1;
opt.NoCheck=0;
opt.Overwrite=1;
origPath=pwd;
cd(pathToData);
if exist('controlRoi+orig.BRIK','file')
    delete controlRoi+orig.BRIK;
    delete controlRoi+orig.HEAD;
end
[err, ErrMessage, Info] = WriteBrik (roiMask, InfoWrite,opt);
cd(origPath);


