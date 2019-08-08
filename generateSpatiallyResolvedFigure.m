% 10.01.18 save t-stats to a BRIK
% 09.22.18 rerunning on p=3 and p=5
% 09.12.18 add computation of AR coefficients
% 09.13.18 add computation of AR residual time series

clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%% load in data
ECHOSTR='3';
FWHMSTR='8';
MAKE_BRAIN_MASK=0;
ALPHA=0.05;
nTRs=645;
NLAGS=5; %p=NLAGS-1 
CONSTANT=0;
LASERONSETTR=215;
LASEROFFSETTR=430;
pathToData='../data/SAVG/NII/';
voxelDim=2.5; % for clustering
minRad=voxelDim;
minClusterVolume=40*voxelDim*voxelDim*voxelDim; % for clustering

% synthesize data filename
optionStr=[ECHOSTR '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2'];
dataFilename=['s' FWHMSTR 'muwBoldEcho' optionStr '+tlrc'];

brainMaskFilename='resampled_brain_mask+tlrc';
[~, smuBold, info, ~] = BrikLoad (fullfile(pathToData,dataFilename));
[~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
brainMask=logical(brainMask);

%% make brain mask FIXTHIS
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
if CONSTANT
    coeffsPreStim=zeros(NLAGS,nVoxels,2); % intercept
    coeffsPrePost=zeros(NLAGS,nVoxels,2);
else
    coeffsPreStim=zeros(NLAGS-1,nVoxels,2); % no intercept
    coeffsPrePost=zeros(NLAGS-1,nVoxels,2);
end
for v=1:nVoxels
    v
    y=boldTs(:,v);
    X=tplitz(y,NLAGS);
    if CONSTANT
        X=X(:,[end 2:end-1]); % include constant term
    else
        X=X(:,2:end-1); % don't include constant term
    end
    
    %% pre versus stim
    [h1,p1,tstat1] = chowtest(X(1:LASEROFFSETTR,:),y(1:LASEROFFSETTR),LASERONSETTR);
    X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
    X2=X(LASERONSETTR+1:LASEROFFSETTR,:); y2=y(LASERONSETTR+1:LASEROFFSETTR);
    coeffsPreStim(:,v,1)=X1\y1;
    coeffsPreStim(:,v,2)=X2\y2;
    
    %% pre versus post
    [h2,p2,tstat2] = chowtest(cat(1,X(1:LASERONSETTR,:),X(LASEROFFSETTR+1:end,:)),cat(1,y(1:LASERONSETTR),y(LASEROFFSETTR+1:end)),LASERONSETTR);
    X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
    X2=X(LASEROFFSETTR+1:end,:); y2=y(LASEROFFSETTR+1:end);
    coeffsPrePost(:,v,1)=X1\y1;
    coeffsPrePost(:,v,2)=X2\y2;
    
    pvals(v,1)=p1; pvals(v,2)=p2;
    tstats(v,1)=tstat1;  tstats(v,2)=tstat2; 
    
end
[pfdr,isSig]=fdr(pvals,ALPHA);

%% make significance volume
% pre vs during
prefix=['isSigBefDur_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc
view='tlrc';
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
prefix=['isSigBefAft_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc

img=ts2vol(isSig(:,2).',brainMask);
[err,ErrMessage,Info]=WriteBrikWrap(pathToData,img,info,prefix,view);
origPath=pwd;
cd(pathToData);
str=['!3dmerge -1clust ' num2str(minRad) ' ' num2str(minClusterVolume) ' -1thresh 0.999 -prefix '  'thresh_' prefix  ' ' prefix '+tlrc'];
eval(str);
cd(origPath);

%% tstat volumes
% before vs during
prefix=['tstatsBefDur_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc
view='tlrc';
img=ts2vol(tstats(:,1).',brainMask);
info.TypeName='float';
info.TypeBytes=4;
info.BRICK_TYPES=3;
[err,ErrMessage,Info]=WriteBrikWrap(pathToData,img,info,prefix,view);

% before vs after
prefix=['tstatsBefAft_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc
view='tlrc';
img=ts2vol(tstats(:,2).',brainMask);
[err,ErrMessage,Info]=WriteBrikWrap(pathToData,img,info,prefix,view);


%%
% resample significance images onto anat grid
cd(pathToData)

threshprefix=['thresh_isSigBefDur_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc
str=['!3dresample -master TT_N27+tlrc -prefix resampled_' threshprefix ' -input ' threshprefix '+tlrc'];
eval(str);

threshprefix=['thresh_isSigBefAft_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc
str=['!3dresample -master TT_N27+tlrc -prefix resampled_' threshprefix ' -input ' threshprefix '+tlrc'];
eval(str);

% resample tstat images onto anat grid
threshprefix=['tstatsBefDur_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc
str=['!3dresample -master TT_N27+tlrc -prefix resampled_' threshprefix ' -input ' threshprefix '+tlrc'];
eval(str);

threshprefix=['tstatsBefAft_NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '_' dataFilename(1:end-5)]; % remove +tlrc
str=['!3dresample -master TT_N27+tlrc -prefix resampled_' threshprefix ' -input ' threshprefix '+tlrc'];
eval(str);


cd(origPath);

%%
% for the sake of analyzing AR spectra, save those
save(['../data/SAVG/NII/arCoeffsResidsEcho' ECHOSTR '-NLAGS' num2str(NLAGS) '_CONSTANT' num2str(CONSTANT) '.mat'],'coeffsPreStim','coeffsPrePost');
