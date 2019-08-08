clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

pathToData='../data/SAVG/NII/';
voxelDim=2.5;
minRad=voxelDim;
minClusterVolume=40*voxelDim*voxelDim*voxelDim; % for clustering

inputPrefix={'isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2'; ...
    'isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho2-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';
    'isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho3-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';
    'isSigBefAft_NLAGS5_CONSTANT0_s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';
    'oc_isSigBefDur_NLAGS5_CONSTANT0_s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';
    'oc_isSigBefAft_NLAGS5_CONSTANT0_s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';}
outputPrefix={'mask_isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2'; ...
    'mask_isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho2-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';
    'mask_isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho3-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';
    'mask_isSigBefAft_NLAGS5_CONSTANT0_s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2'
    'mask_oc_isSigBefDur_NLAGS5_CONSTANT0_s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';
    'mask_oc_isSigBefAft_NLAGS5_CONSTANT0_s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2';};

nFiles=numel(inputPrefix);

for f=1:nFiles
    origPath=pwd;
    cd(pathToData);
    str=['!3dmerge -1clust_order ' num2str(minRad) ' ' num2str(minClusterVolume) ' -1thresh 0.999 -prefix '  outputPrefix{f}  ' ' inputPrefix{f} '+tlrc'];
    eval(str);
    cd(origPath);
end

%%
% resample onto anatomy
cd(pathToData)
for f=1:nFiles
    str=['!3dresample -master TT_N27+tlrc -prefix ' ['resampled_' outputPrefix{f}] ' -input ' outputPrefix{f} '+tlrc'];
    eval(str);
end
cd(origPath)