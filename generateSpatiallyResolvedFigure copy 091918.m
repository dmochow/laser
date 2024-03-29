% 09.12.18 add computation of AR coefficients
% 09.13.18 add computation of AR residual time series

clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%% load in data
MAKE_BRAIN_MASK=0;
ALPHA=0.05;
FWHMSTR='8';
ECHOSTR='3';
nTRs=645;
NLAGS=5;
LASERONSETTR=215;
LASEROFFSETTR=430;
WHITE_MATTER=1; CSF=0; Kwm=3; Kcsf=1; CENTERBOLDS=1;
pathToData='../data/SAVG/NII/';
dataFilename=['s' FWHMSTR 'muwBoldEcho' ECHOSTR '-WhiteMatter' ... 
    num2str(WHITE_MATTER) '-KWHM' num2str(Kwm) '-CSF' num2str(CSF) ...
    '-KCSF' num2str(Kcsf) '-Center' num2str(CENTERBOLDS) '+tlrc'];

brainMaskFilename='resampled_brain_mask+tlrc';
[~, smuBold, info, ~] = BrikLoad (fullfile(pathToData,dataFilename));
[~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
brainMask=logical(brainMask);
voxelDim=2.5; % for clustering
minRad=voxelDim;
minClusterVolume=40*voxelDim*voxelDim*voxelDim; % for clustering

%% make brain mask
if MAKE_BRAIN_MASK
    !export PATH=$PATH:/Users/jacekdmochowski/abin
    PATH = getenv('PATH');
    setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
    origPath=pwd;
    cd(pathToData);
    str='!3dAutomask -prefix brain_mask TT_N27+tlrc';
    eval(str);
    str2=['!3dresample -master ' dataFilename ' -prefix resampled_brain_mask -input brain_mask+tlrc'];
    eval(str2);
    cd(origPath);
end

%%
boldTs=vol2ts(smuBold,brainMask);
nVoxels=size(boldTs,2);
pvals=zeros(nVoxels,2);
tstats=zeros(nVoxels,2);
%coeffsPreStim=zeros(NLAGS,nVoxels,2); % intercept
%coeffsPrePost=zeros(NLAGS,nVoxels,2);
coeffsPreStim=zeros(NLAGS-1,nVoxels,2); % no intercept
coeffsPrePost=zeros(NLAGS-1,nVoxels,2);
rsdPreStim=zeros(nTRs,nVoxels,2);
rsdPrePost=zeros(nTRs,nVoxels,2);
for v=1:nVoxels
    v
    y=boldTs(:,v);
    X=tplitz(y,NLAGS); 
    %X=X(:,[end 2:end-1]); % include constant term
    X=X(:,2:end-1); % don't include constant term
    
    %% pre versus stim
    [h1,p1,tstat1] = chowtest(X(1:LASEROFFSETTR,:),y(1:LASEROFFSETTR),LASERONSETTR);
    %[h1,p1,tstat1,cValue,SS] = mychowtest(X(1:LASEROFFSETTR,:),y(1:LASEROFFSETTR),LASERONSETTR);
    X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
    X2=X(LASERONSETTR+1:LASEROFFSETTR,:); y2=y(LASERONSETTR+1:LASEROFFSETTR);
    coeffsPreStim(:,v,1)=X1\y1; 
    coeffsPreStim(:,v,2)=X2\y2; 
    rsdPreStim(:,v,1)=y-X*coeffsPreStim(:,v,1);
    rsdPreStim(:,v,2)=y-X*coeffsPreStim(:,v,2);
    
    %% pre versus post
    [h2,p2,tstat2] = chowtest(cat(1,X(1:LASERONSETTR,:),X(LASEROFFSETTR+1:end,:)),cat(1,y(1:LASERONSETTR),y(LASEROFFSETTR+1:end)),LASERONSETTR);
    pvals(v,1)=p1; pvals(v,2)=p2;
    tstats(v,1)=tstat1; tstats(v,2)=tstat2;
    X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
    X2=X(LASEROFFSETTR+1:end,:); y2=y(LASEROFFSETTR+1:end);
    coeffsPrePost(:,v,1)=X1\y1; 
    coeffsPrePost(:,v,2)=X2\y2; 
    rsdPrePost(:,v,1)=y-X*coeffsPrePost(:,v,1);
    rsdPrePost(:,v,2)=y-X*coeffsPrePost(:,v,2);
end
[pfdr,isSig]=fdr(pvals,ALPHA);

%% make significance volume
% pre vs during
prefix=['isSigBefDur_FWHM' FWHMSTR '_ECHO' ECHOSTR '-WhiteMatter' ... 
    num2str(WHITE_MATTER) '-KWM' num2str(Kwm) '-CSF' num2str(CSF) ...
    '-KCSF' num2str(Kcsf) '-Center' num2str(CENTERBOLDS)]; view='tlrc'; 
img=ts2vol(isSig(:,1).',brainMask);
[err,ErrMessage,Info]=WriteBrikWrap(pathToData,img,info,prefix,view);
% cluster and prune
origPath=pwd;
cd(pathToData);
str=['!3dmerge -1clust ' num2str(minRad) ' ' num2str(minClusterVolume) ' -1thresh 0.999 -prefix '  'thresh_' prefix  ' ' prefix '+tlrc'];
eval(str);
cd(origPath);
%%
% pre vs after
prefix=['isSigBefAft_FWHM' FWHMSTR '_ECHO' ECHOSTR '-WhiteMatter' ... 
    num2str(WHITE_MATTER) '-KWM' num2str(Kwm) '-CSF' num2str(CSF) ...
    '-KCSF' num2str(Kcsf) '-Center' num2str(CENTERBOLDS)]; view='tlrc'; 
img=ts2vol(isSig(:,2).',brainMask);
[err,ErrMessage,Info]=WriteBrikWrap(pathToData,img,info,prefix,view);
origPath=pwd;
cd(pathToData);
str=['!3dmerge -1clust ' num2str(minRad) ' ' num2str(minClusterVolume) ' -1thresh 0.999 -prefix '  'thresh_' prefix  ' ' prefix '+tlrc'];
eval(str);
cd(origPath);
%%
% resample significance images onto anat grid
cd(pathToData)

% before vs during
threshprefix=['thresh_isSigBefDur_FWHM' FWHMSTR '_ECHO' ECHOSTR '-WhiteMatter' ... 
    num2str(WHITE_MATTER) '-KWM' num2str(Kwm) '-CSF' num2str(CSF) ...
    '-KCSF' num2str(Kcsf) '-Center' num2str(CENTERBOLDS)];
str=['!3dresample -master TT_N27+tlrc -prefix resampled_' threshprefix ' -input ' threshprefix '+tlrc'];
eval(str);

% before vs after
threshprefix=['thresh_isSigBefAft_FWHM' FWHMSTR '_ECHO' ECHOSTR '-WhiteMatter' ... 
    num2str(WHITE_MATTER) '-KWM' num2str(Kwm) '-CSF' num2str(CSF) ...
    '-KCSF' num2str(Kcsf) '-Center' num2str(CENTERBOLDS)];
str=['!3dresample -master TT_N27+tlrc -prefix resampled_' threshprefix ' -input ' threshprefix '+tlrc'];
eval(str);

cd(origPath);

%%
% for the sake of correlating absorption to effect, write out tstats to
%mat volumes
vstat1=ts2vol(tstats(:,1).',brainMask);
vstat2=ts2vol(tstats(:,2).',brainMask);
vpvals1=ts2vol(pvals(:,1).',brainMask);
vpvals2=ts2vol(pvals(:,2).',brainMask);
save('../data/SAVG/NII/tstatVolume.mat','vstat1','vstat2','vpvals1','vpvals2');

%%
% for the sake of analyzing AR spectra, save those
save(['../data/SAVG/NII/arCoeffsResidsEcho' ECHOSTR '.mat'],'coeffsPreStim','coeffsPrePost',...
    'rsdPreStim','rsdPrePost');

