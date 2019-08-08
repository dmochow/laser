clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin

% update 02/13/18
% work in cylindrical coordinates
% left-right along laser face is the reference direction rho
% z is perpendicular to that, namely depth into the head
%
% 02/13/18 PM: rotate axis to align laser face with reference plane

% params go here
% NB
% dim 1: sagittal (L->R)
% dim 2: coronal (P->A) (afni display is A->P)
% dim 3: axial (I->S)

subjStr='S10';
ulf=[141 23 181]+[1 1 1]; 
urf=[158 28 177]+[1 1 1]; 
urb=[165 9 180]+[1 1 1]; 
lrf=[155 18 137]+[1 1 1];

% subjStr='S11';
% ulf=[127 48 197]+[1 1 1]; 
% urf=[148 55 193]+[1 1 1]; 
% urb=[161 25 213]+[1 1 1]; 
% lrf=[149 30 158]+[1 1 1];

% try to compute brain mask in afni
% THIS DOESN'T WORK DUE TO MATLAB NOT SEEING THE BASH PATH
% pathToData=['../data/' subjStr '/NII/'];
% origPath=pwd;
% cd(pathToData);
% str='!3dSkullStrip anat.nii';
% eval(str);
% cd(origPath);

% processing begins here
pathToData=['../data/' subjStr '/NII/'];
niiFilename=['../data/' subjStr '/NII/anat.nii'];
brainMaskFilename=['../data/' subjStr '/NII/brain_mask+orig'];
outNiifilename=['../data/' subjStr '/NII/roiMask.nii'];

laserOrigin = getLaserOrigin(ulf,urf,urb); % origin in afni space, where second dimension is A->P

%%
[err, mask, Info, ErrMessage] = BrikLoad (brainMaskFilename);
mridim=size(mask);
[x,y,z]=ndgrid(0:mridim(1)-1,0:mridim(2)-1,0:mridim(3)-1);
LOx=laserOrigin(1)*ones(size(x));
LOy=laserOrigin(2)*ones(size(y));
LOz=laserOrigin(3)*ones(size(z));

%% 
% rotation matrices
laserOrigin(2)=256-laserOrigin(2);
theta=atan(laserOrigin(1)/laserOrigin(2));
phi=asin(laserOrigin(3)/norm(laserOrigin));
%theta*180/pi
%phi*180/pi

%%
Rz=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
Rx=[1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi) ];
R=Rz*Rx;

% translation
xp=x-laserOrigin(1);
yp=y-laserOrigin(2);
zp=z-laserOrigin(3);

% rotation
xp=R(1,1)*xp+R(1,2)*yp+R(1,3)*zp;
yp=R(2,1)*xp+R(2,2)*yp+R(2,3)*zp;
zp=R(3,1)*xp+R(3,2)*yp+R(3,3)*zp;

% coordinates of laser center in new space
xp(round(laserOrigin(1)),round(laserOrigin(2)),round(laserOrigin(3)))
yp(round(laserOrigin(1)),round(laserOrigin(2)),round(laserOrigin(3)))
zp(round(laserOrigin(1)),round(laserOrigin(2)),round(laserOrigin(3)))

% radius=sqrt((xx-LOx).^2+(yy-LOy).^2+(zz-LOz).^2);
% azimuth=atan2d( (xx-LOx),(yy-LOy) );
% elevation=acosd((zz-LOz)./radius)-90;
% 
% % compute rho and z
% z=(xx-LOx).*crossprod(1)+(yy-LOy).*crossprod(2)+(zz-LOz).*crossprod(3); % define left right
% rho=(xx-LOx).*lrvec(1)+(yy-LOy).*lrvec(2)+(zz-LOz).*lrvec(3); % define front back

roiMask=zeros(size(mask));
% roiMask=tmp;
%roiMask(radius<30 & abs(azimuth)<45 & abs(elevation)<45 )=1;
roiMask(abs(zp)<30)=1;

%% write out radial mask
InfoWrite=Info;
opt.Scale=0;
opt.Prefix='roi';
opt.View='orig';
opt.Verbose=1;
opt.AppendHistory=1;
opt.NoCheck=0;
opt.Overwrite=1;
roiMask=roiMask(:,end:-1:1,:);  % apparently needed to avoid flipping PA dim
origPath=pwd;
cd(pathToData);
[err, ErrMessage, Info] = WriteBrik (roiMask, InfoWrite,opt);
cd(origPath);


% %% 
% % run the search
% dr=0.1; % radial increment
% currentPoint=laserOrigin;
% maxIter=1e6;
% 
% crossprod=cross((ulf-urf),(ulf-lrf));
% crossprod=crossprod/norm(crossprod);
% 
% % if pointing to the front (anterior), flip it
% if crossprod(2)<0
%     crossprod(2)=-crossprod(2);
% end
% 
% for i=1:maxIter
%     currentPoint=currentPoint+dr*crossprod;
%     xq=round(currentPoint); % query point which must be discrete!
%     if mask(xq(1)+1,mridim(2)-xq(2),xq(3)+1)
%         break
%     else
%         continue
%     end
% end
% 
% nearestBrainPoint=xq;
% 
% %% now figure out the parallel planes on each side of the nearestBrainPoint
% lrvec=urf-ulf;
% lrvec=lrvec/norm(lrvec);
% if lrvec(1)>0
%     leftNearestBrainPoint=nearestBrainPoint-dLeft*lrvec; % move left 
%     rightNearestBrainPoint=nearestBrainPoint+dLeft*lrvec; 
% else
%     leftNearestBrainPoint=nearestBrainPoint+dLeft*lrvec; 
%     rightNearestBrainPoint=nearestBrainPoint-dLeft*lrvec; 
% end
% 
% nearestBrainPoint
% leftNearestBrainPoint
% rightNearestBrainPoint


%% parallel planes on top/bottom of nearestBrainPoint

%%
% for debugging the coordinate transformation between Afni and Matlab via
% BrikLoad
% figure; 
% subplot(221);
% imagesc(squeeze(mask(x+1,:,:)));
% subplot(222);
% imagesc(squeeze(mask(:,mridim(2)-y,:)));  % don't ask, it works
% subplot(223);
% imagesc(squeeze(mask(:,:,z+1)));





