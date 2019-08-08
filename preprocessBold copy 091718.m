
% 08/28/18b: remove projection of absorption
% add talairaching and resampling of laser origin
% 08/28/18
% added projection of absorption + talairach + resampling on bold
% 07/17/18
% modularizing preprocessBold:
% brain mask generation
% bold-epi+brain mask registration + resampling
% roi generation + registration + resampling
% 07/25/18
% subjects without markers handled by not doing any ROI generation or
% processing
% the chain is now:
% (1) preprocessBold
% (2) makeAverageROI
% (3) postprocessBold (tentatively)
% UPDATED (08/28/18) the chain is now:
% (1) preprocessBold
% (2) postprocessBold 
% (3) alignAbsorption


clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

PROCESS_ONLY_ROI=1;
%subjStrs={'S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17',...
 %    'S18','S19','S20','S21','S22'};
subjStrs={'S11'};
nSubjects=numel(subjStrs);
origPath=pwd;
%roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
%controlRoiPrefix=['control_roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
%Aprefix=['A-' date]; % the prefix for the absorption data
laserOriginPrefix='laserOriginV';

for s=1:nSubjects
    subjStr=subjStrs{s};
    [urf,ulf,urb,lrf,midsag]=getMarkerCoords(subjStr);
    pathToData=['../data/' subjStr '/NII/'];
    anatFilename=fullfile(pathToData,'anat+orig'); %['../data/' subjStr '/NII/anat+orig'];
    anatWithSkullFilename=fullfile(pathToData,'anatWithSkull+orig'); %['../data/' subjStr '/NII/anatWithSkull+orig'];
    brainMaskFilename=fullfile(pathToData,'brain_mask+orig'); %['../data/' subjStr '/NII/brain_mask+orig'];
    
    if ~PROCESS_ONLY_ROI % only do this if we want to generate everything from scratch
        
        %% delete all but the original anat.nii, bold_ex.nii files
        filenames=dir(pathToData);
        for f=1:numel(filenames)
            if ~ (strcmp(filenames(f).name,'anat.nii') || strcmp(filenames(f).name,'bold_e1.nii') ...
                    || strcmp(filenames(f).name,'bold_e2.nii') || strcmp(filenames(f).name,'bold_e3.nii'))
                if isfile(fullfile(pathToData,filenames(f).name))
                    eval(['delete ' fullfile(pathToData,filenames(f).name)]);
                end
            end
        end
        %%
        
        %% Part I: brain mask generation
        cd(pathToData);
        if ~exist('brain_mask+orig.HEAD','file')
            str1='!3dcopy anat.nii anat';
            str1b='!3dcopy anat.nii anatWithSkull';
            str2='!3dSkullStrip -input anat+orig';
            str3='!3dAutomask -prefix brain_mask skull_strip_out+orig';
            str3b='!gunzip brain_mask+orig.BRIK.gz';
            eval(str1);
            eval(str1b);
            eval(str2);
            eval(str3);
            eval(str3b);
        end
        cd(origPath);
        
        %% Part II: bold-anatomy registration, atlasing, and resampling
        cd(pathToData);
        str4='!3dcopy bold_e1.nii bold_e1'; eval(str4);
        str5='!3dcopy bold_e2.nii bold_e2'; eval(str5);
        str6='!3dcopy bold_e3.nii bold_e3'; eval(str6);
        str7='!3ddespike -NEW -localedit -prefix dsbold_e1 bold_e1+orig'; eval(str7);
        str8='!3ddespike -NEW -localedit -prefix dsbold_e2 bold_e2+orig'; eval(str8);
        str9='!3ddespike -NEW -localedit -prefix dsbold_e3 bold_e3+orig'; eval(str9);
        % align the anatomical to the talairach atlas
        str10='!@auto_tlrc -base TT_N27+tlrc -input anat+orig'; eval(str10);
        % align the brain mask to the talairach atlas
        str10b='!@auto_tlrc -apar anat+tlrc -input brain_mask+orig'; eval(str10b);
        % now align the epis to the anatomy
        str12='!align_epi_anat.py -anat anat+orig -epi dsbold_e2+orig -epi_base 5 -child_epi dsbold_e1+orig dsbold_e3+orig -epi2anat -suffix _al -tlrc_apar anat+tlrc';
        eval(str12);
        % resample the brain mask on the bold grid
        str14b='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_brain_mask -input brain_mask+tlrc'; eval(str14b);
        cd(origPath);
        
    end
    
    %% Part III: roi generation, atlasing, and resampling
    deleteRoiFiles(pathToData);
 
    if ~(isempty(urf)||isempty(ulf)||isempty(lrf)||isempty(urb))
    
    laserOrigin = getLaserOrigin(ulf,urf,urb);
    laserOrigin = fineTuneOrigin(anatWithSkullFilename,laserOrigin,ulf,urf,lrf);
    save(fullfile(pathToData,'laserOrigin.mat'),'laserOrigin');
    
    % make laser origin BRIK here
    [~, tmp, Info, ~] = BrikLoad (anatFilename); % only to get the 'Info' field  and size of MRI 
    laserOriginV=zeros(size(tmp));
    laserOriginV(round(laserOrigin(1)),round(laserOrigin(2)),round(laserOrigin(3)))=1;
    WriteBrikWrap(pathToData,laserOriginV,Info,laserOriginPrefix,'orig');

    
%     Amri = projectAbsorption(brainMaskFilename,laserOrigin,urf,ulf,urb,lrf);
%     %roiMask = projectLaser(brainMaskFilename,laserOrigin,urf,ulf,urb,lrf,radMask,depthMask);  % project the real laser into the mri
%     [~, ~, Info, ~] = BrikLoad (anatFilename); % only to get the 'Info' field  
%     Info.TypeName='float';
%     Info.TypeBytes=4;
%     Info.BRICK_TYPES=3;
%     WriteBrikWrap(pathToData,Amri,Info,Aprefix,'orig');

    % now resample laser origin mri
    cd(pathToData);
    str11=['!@auto_tlrc -apar anat+tlrc -input ' laserOriginPrefix '+orig'] ; eval(str11); % register A to anatomy
    str11a=['!gunzip ' laserOriginPrefix '+tlrc.BRIK.gz']; eval(str11a); % weird that this is needed
    str13=['!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_' laserOriginPrefix ' -input ' laserOriginPrefix '+tlrc']; eval(str13);  % resample on bold grid
    
    cd(origPath);
    
    end
end