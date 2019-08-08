% 08/29/18: project absorption onto standard head
clear all; close all; clc

pathToData='../data/SAVG/NII';
brainMaskFilename=fullfile(pathToData,'brain_mask+tlrc');
subjStrs={'S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17',...
    'S18','S19','S20','S21','S22'}; % all subjects with markers
nSubjects=numel(subjStrs);



%% need laserOrigin, urf,... so...
% voxel space to mri space
% mri space to talairach mri space
% talairach mri space to voxel space
XYZOUT=zeros(3,nSubjects);
URFOUT=zeros(3,nSubjects);
ULFOUT=zeros(3,nSubjects);
URBOUT=zeros(3,nSubjects);
LRFOUT=zeros(3,nSubjects);

for s=1:nSubjects
    subjStr=subjStrs{s};
    % get corners of headhear
    [urf,ulf,urb,lrf,midsag]=getMarkerCoords(subjStr);
    pathToData=['../data/' subjStr '/NII/'];
    
    % get laser origin
    load(fullfile(pathToData,'laserOrigin.mat'),'laserOrigin');
    Xin=[laserOrigin';1];
    Xin=round(Xin);
    
    % get talairach transform matrix Xat
    xatFilename=fullfile(pathToData,'anat.Xat.1D');
    Xat = dlmread(xatFilename);
    Xat=cat(1,Xat,[0 0 0 1]);
    
    % get headers
    [~,~,InfoOrig,~]=BrikLoad(fullfile(pathToData,'anat+orig.BRIK'));
    [~,~,InfoTlrc,~]=BrikLoad(fullfile(pathToData,'anat+tlrc.BRIK'));
    
    % transform
    [xyzOut,voxOut] = orig2tal(Xin,Xat,InfoOrig,InfoTlrc);
    [urfOut,urfOutVox] = orig2tal(urf(:),Xat,InfoOrig,InfoTlrc);
    [ulfOut,ulfOutVox] = orig2tal(ulf(:),Xat,InfoOrig,InfoTlrc);
    [urbOut,urbOutVox] = orig2tal(urb(:),Xat,InfoOrig,InfoTlrc);
    [lrfOut,lrfOutVox] = orig2tal(lrf(:),Xat,InfoOrig,InfoTlrc);
    
    XYZOUT(:,s)=xyzOut;
    URFOUT(:,s)=urfOut;
    ULFOUT(:,s)=ulfOut;
    URBOUT(:,s)=urbOut;
    LRFOUT(:,s)=lrfOut;
    
    
    
    
end

muXYZ=mean(XYZOUT,2);
muURF=mean(URFOUT,2);
muULF=mean(ULFOUT,2);
muURB=mean(URBOUT,2);
muLRF=mean(LRFOUT,2);

stLaserOrigin = xyz2vox(muXYZ,InfoTlrc);
stUrf= xyz2vox(muURF,InfoTlrc);
stUlf= xyz2vox(muULF,InfoTlrc);
stUrb= xyz2vox(muURB,InfoTlrc);
stLrf= xyz2vox(muLRF,InfoTlrc);

% finally project
Amri = projectAbsorption(brainMaskFilename,stLaserOrigin,stUrf,stUlf,stUrb,stLrf);
[~, ~, Info, ~] = BrikLoad (brainMaskFilename); % only to get the 'Info' field
Info.TypeName='float';
Info.TypeBytes=4;
Info.BRICK_TYPES=3;
WriteBrikWrap('../data/SAVG/NII',Amri,Info,'Amu','tlrc');

%% ... and resample (hopefully the last bit of code for this project)
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
origPath=pwd;
cd('../data/SAVG/NII');
str=['!3dresample -master resampled_brain_mask+tlrc -prefix resampled_Amu -input Amu+tlrc']; 
eval(str); 
cd(origPath);

%% FROM SCRIPT
%  % from voxel orig to XYZ orig
%     Xout=zeros(3,1);
%     Xout(1)=InfoOrig.ORIGIN(1)+Xin(1)*InfoOrig.DELTA(1);
%     Xout(2)=InfoOrig.ORIGIN(2)+Xin(2)*InfoOrig.DELTA(2);
%     Xout(3)=InfoOrig.ORIGIN(3)+Xin(3)*InfoOrig.DELTA(3);
%     Xout=[Xout(:);1];
%     
%     % from XYZ orig to XYZ talairach
%     Xout2=inv(Xat)*Xout;
%     Xout2=Xout2(1:3);
%     Xout2=Xout2(:);
%     
%     
%     % from XYZ talairach to voxel talairach
%     Xout3=zeros(3,1);
%     Xout3(1)= (Xout2(1)-InfoTlrc.ORIGIN(1))/InfoTlrc.DELTA(1);
%     Xout3(2)= (Xout2(2)-InfoTlrc.ORIGIN(2))/InfoTlrc.DELTA(2);
%     Xout3(3)= (Xout2(3)-InfoTlrc.ORIGIN(3))/InfoTlrc.DELTA(3);

%% SCRAPHEAP
%
% %Xin=[109,205,153]; % most anterior point of brain
% 
% [~,anat,InfoOrig,~]=BrikLoad('../data/S11/NII/anat+orig.BRIK');
% [~,anatTlrc,InfoTlrc,~]=BrikLoad('../data/S11/NII/anat+tlrc.BRIK');
% Xat=[1.03085,0.00213217,-0.013939,-3.20999;
%     0.0420256,1.00926,0.012311,-36.0469;
%     -0.0148486,-0.136752,1.01311,-31.2451;
%     0,0,0,1];
% Xorig=reshape(InfoOrig.IJK_TO_DICOM_REAL,4,3)';
% Xorig=cat(1,Xorig,[0 0 0 1]);
% Xtlrc=reshape(InfoTlrc.IJK_TO_DICOM_REAL,4,3)';
% Xtlrc=cat(1,Xtlrc,[0 0 0 1]);
% 
% % XYZ orig
% Xout(1)=InfoOrig.ORIGIN(1)+Xin(1)*InfoOrig.DELTA(1);
% Xout(2)=InfoOrig.ORIGIN(2)+Xin(2)*InfoOrig.DELTA(2);
% Xout(3)=InfoOrig.ORIGIN(3)+Xin(3)*InfoOrig.DELTA(3);
% Xout=[Xout(:);1];
% 
% 
% % from XYZ orig to XYZ talairach
% Xout2=inv(Xat)*Xout;
% Xout2=Xout2(1:3);
% Xout2=Xout2(:);
% 
% Xout3(1)= (Xout2(1)-InfoTlrc.ORIGIN(1))/InfoTlrc.DELTA(1);
% Xout3(2)= (Xout2(2)-InfoTlrc.ORIGIN(2))/InfoTlrc.DELTA(2);
% Xout3(3)= (Xout2(3)-InfoTlrc.ORIGIN(3))/InfoTlrc.DELTA(3);
% 
% 
% 
% % % origin of orig to origin of tal
% % originOrig=cat(1,InfoOrig.ORIGIN(:),1);
% % originTal=cat(1,InfoTlrc.ORIGIN(:),1);
% % % to go from orig to tlrc: tlrc=inv(Xat)*orig
% %
% %
% % y=Xat*x;
% % z=y; z(2)=255-y(2);
% % [x y z];
