% 05/30/18
% nuisance regression and smoothing taken out as it is now handled in
% postprocessBold.m

% we need the following files present in the pathToData folder:
% anat.nii
% bold_e1.nii
% bold_e2.nii
% bold_e3.nii

clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

subjStrs={'S22'};
ONLY_PROCESS_ROI=1;
%subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'};
nSubjects=numel(subjStrs);

for s=1:nSubjects
    subjStr=subjStrs{s};
    pathToData=['../data/' subjStr '/NII/'];
    
    switch subjStr
        case 'S02'
        case 'S03'
        case 'S04'
            urf=[158 255-26 164]+[1 1 1];
            ulf=[123 255-24 166]+[1 1 1];
            urb=[164 255-8 164]+[1 1 1];
            lrf=[156 255-28 151]+[1 1 1];
            midsag=113;
        case 'S05'
            urf=[150 255-22 159]+[1 1 1];
            ulf=[116 255-20 161]+[1 1 1];
            urb=[156 255-1 164]+[1 1 1];
            lrf=[150 255-22 151]+[1 1 1];
            midsag=113;
        case 'S06'
            urf=[161 255-20 170]+[1 1 1];
            ulf=[128 255-17 170]+[1 1 1];
            urb=[169 255-2 170]+[1 1 1];
            lrf=[161 255-19 162]+[1 1 1];
            midsag=112;
        case 'S07'
            urf=[158 255-27 165]+[1 1 1];
            ulf=[126 255-20 171]+[1 1 1];
            urb=[165 255-6 167]+[1 1 1];
            lrf=[156 255-27 156]+[1 1 1];
            midsag=113;
        case 'S09'
            ulf=[126 255-43 179]+[1 1 1];
            urf=[147 255-49 179]+[1 1 1];
            urb=[154 255-35 192]+[1 1 1];
            lrf=[138 255-17 149]+[1 1 1];
            midsag=113;
        case 'S10'
            ulf=[141 255-23 181]+[1 1 1];
            urf=[158 255-28 177]+[1 1 1];
            urb=[165 255-9 180]+[1 1 1];
            lrf=[155 255-18 137]+[1 1 1];
            midsag=113;
        case 'S11'
            ulf=[127 255-48 197]+[1 1 1];
            urf=[148 255-55 193]+[1 1 1];
            urb=[161 255-25 213]+[1 1 1];
            lrf=[149 255-30 158]+[1 1 1];
            midsag=111;
        case 'S12'
            urf=[151 255-47 198]+[1 1 1]; %
            ulf=[129 255-41 202]+[1 1 1]; %
            lrf=[151 255-34 159]+[1 1 1]; %
            urb=[152 255-32 203]+[1 1 1]; %
            midsag=112;
        case 'S13'
            urf=[146 255-31 182]+[1 1 1]; %
            ulf=[123 255-28 186]+[1 1 1]; %
            urb=[153 255-3 192]+[1 1 1];
            lrf=[144 255-19 143]+[1 1 1]; %
            midsag=112;
        case 'S14'
            urf=[141 255-44 192]+[1 1 1];
            ulf=[122 255-39 192]+[1 1 1];
            urb=[148 255-23 197]+[1 1 1];
            lrf=[150 255-40 152]+[1 1 1];
            midsag=110;
        case 'S15'
            urf=[155 255-45 178]+[1 1 1];
            ulf=[137 255-35 178]+[1 1 1];
            urb=[171 255-14 189]+[1 1 1];
            lrf=[165 255-36 139]+[1 1 1];
            midsag=113;
        case 'S16'
            urf=[152 255-44 159]+[1 1 1]; %
            ulf=[132 255-36 161]+[1 1 1]; %
            lrf=[154 255-36 118]+[1 1 1]; %
            urb=[169 255-9 167]+[1 1 1]; %
            midsag=112;
        case 'S17'
            urf=[149 255-32 155]+[1 1 1]; %
            ulf=[130 255-26 155]+[1 1 1]; %
            lrf=[152 255-34 112]+[1 1 1]; %
            urb=[157 255-4 153]+[1 1 1]; %
            midsag=112;
        case 'S18'
            urf=[144 255-36 198]+[1 1 1]; %
            ulf=[126 255-31 197]+[1 1 1]; %
            lrf=[150 255-22 152]+[1 1 1]; %
            urb=[153 255-9 202]+[1 1 1]; %
            midsag=112;          
        case 'S19'
            urf=[144 255-48 200]+[1 1 1]; %
            ulf=[125 255-43 200]+[1 1 1]; %
            lrf=[147 255-34 161]+[1 1 1]; %
            urb=[153 255-14 212]+[1 1 1]; %
            midsag=112;
        case 'S20'
            urf=[156 255-44 193]+[1 1 1]; %
            ulf=[139 255-36 193]+[1 1 1]; %
            lrf=[162 255-37 153]+[1 1 1]; %
            urb=[167 255-18 199]+[1 1 1]; %
            midsag=112;          
        case 'S21'
            urf=[151 255-48 203]+[1 1 1]; %
            ulf=[134 255-42 202]+[1 1 1]; %
            lrf=[156 255-29 161]+[1 1 1]; %
            urb=[161 255-25 212]+[1 1 1]; %
            midsag=113;         
        case 'S22'
            urf=[139 255-50 206]+[1 1 1]; %
            ulf=[125 255-48 206]+[1 1 1]; %
            lrf=[144 255-34 168]+[1 1 1]; %
            urb=[145 255-26 217]+[1 1 1]; %
            midsag=113;
    end

    niiFilename=fullfile(pathToData,'anat.nii');
    anatFilename=fullfile(pathToData,'anat+orig');
    brainMaskFilename=fullfile(pathToData,'brain_mask+orig');
    outNiiFilename=fullfile(pathToData,'roiMask.nii');
    outControlNiiFilename=fullfile(pathToData,'controlRoiMask.nii');
    
    % compute brain mask in afni
    origPath=pwd;
    cd(pathToData);
    if ~exist('brain_mask+orig.HEAD','file')
        str1='!3dcopy anat.nii anat';
        str2='!3dSkullStrip -input anat+orig';
        str3='!3dAutomask -prefix brain_mask skull_strip_out+orig';
        eval(str1);
        eval(str2);
        eval(str3);
    end
    cd(origPath);
    
    hasMarkers=0;
    if exist('ulf') && exist('urf') && exist ('urb') && exist ('lrf')
        hasMarkers=1;
        laserOrigin = getLaserOrigin(ulf,urf,urb);
        
        % reflect the laser origin to get the control (contralateral) laser
        % origin
        urfc=urf; ulfc=ulf; urbc=urb; lfrc=lrf;
        urfc(1)=2*midsag-urf(1);
        ulfc(1)=2*midsag-ulf(1);
        urbc(1)=2*midsag-urb(1);
        lrfc(1)=2*midsag-lrf(1);
        %         controlLaserOrigin = getLaserOrigin(ulfc,urfc,urbc);
        
        % fine-tune the position of real laser (project it onto the scalp)
        laserOrigin = fineTuneOrigin(anatFilename,laserOrigin,ulf,urf,lrf);
        save(fullfile(pathToData,'laserOrigin.mat'),'laserOrigin');
        
        % make the control ROI
        controlLaserOrigin=laserOrigin; % AFTER fine-tuning
        controlLaserOrigin(1)=2*midsag-laserOrigin(1);
        save(fullfile(pathToData,'controlLaserOrigin.mat'),'controlLaserOrigin');
        
        % project the real laser into the mri
        roiMask = projectLaser(brainMaskFilename,laserOrigin,urf,ulf,urb,lrf);
        [~, ~, Info, ~] = BrikLoad (anatFilename); % we only do this to get the 'Info' field
        WriteBrikWrap(pathToData,roiMask,Info,'roi','orig');
        
        % project the control laser onto the mri
        controlRoiMask = projectLaser(brainMaskFilename,controlLaserOrigin,urfc,ulfc,urbc,lrfc);
        WriteBrikWrap(pathToData,controlRoiMask,Info,'controlRoi','orig');
        
    end
    
    %%
    % entering what used to be 'runAfniPreProcessing.m'
    if ~ONLY_PROCESS_ROI
        cd(pathToData);
        str4='!3dcopy bold_e1.nii bold_e1';
        eval(str4);
        
        str5='!3dcopy bold_e2.nii bold_e2';
        eval(str5);
        
        str6='!3dcopy bold_e3.nii bold_e3';
        eval(str6);
        
        str7='!3ddespike -NEW -localedit -prefix dsbold_e1 bold_e1+orig';
        eval(str7);
        
        str8='!3ddespike -NEW -localedit -prefix dsbold_e2 bold_e2+orig';
        eval(str8);
        
        str9='!3ddespike -NEW -localedit -prefix dsbold_e3 bold_e3+orig';
        eval(str9);
        
        % align the anatomical to the talairach atlas
        str10='!@auto_tlrc -base TT_N27+tlrc -input anat+orig';
        eval(str10);
        
        % align the brain mask (derived from skull stripped anatomical) with
        % the talairach atlas
        str10b='!@auto_tlrc -apar anat+tlrc -input brain_mask+orig';
        eval(str10b);
        
        % now align the epis to the anatomy, and also telling it that we have an
        % atlased anatomy
        str12='!align_epi_anat.py -anat anat+orig -epi dsbold_e2+orig -epi_base 5 -child_epi dsbold_e1+orig dsbold_e3+orig -epi2anat -suffix _al -tlrc_apar anat+tlrc';
        eval(str12);
        
        % commenting out because we are doing this in postprocessBold.m
        %         str12a='!3dTproject -input dsbold_e1_tlrc_al+tlrc -prefix nudsbold_e1_tlrc_al -ort dsbold_e1_vr_motion.1D';
        %         eval(str12a);
        %
        %         str12b='!3dTproject -input dsbold_e2_tlrc_al+tlrc -prefix nudsbold_e2_tlrc_al -ort dsbold_e2_vr_motion.1D';
        %         eval(str12b);
        %
        %         str12c='!3dTproject -input dsbold_e3_tlrc_al+tlrc -prefix nudsbold_e3_tlrc_al -ort dsbold_e3_vr_motion.1D';
        %         eval(str12c);
        
        
        if hasMarkers
            % align the roi mask to the talairach atlas
            str11='!@auto_tlrc -apar anat+tlrc -input roi+orig'; % this produces roi+tlrc
            eval(str11);
            
            str11b='!@auto_tlrc -apar anat+tlrc -input controlRoi+orig'; % produces controlRoi+tlrc
            eval(str11b);
            
            % resample the ROI on the bold grid
            %             str13='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_roi -input roi+tlrc';
            str13='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_roi -input roi+tlrc';
            eval(str13);
            
            %             str13b='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_control_roi -input controlRoi+tlrc';
            str13b='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_control_roi -input controlRoi+tlrc';
            eval(str13b);
        end
        
        % resample the brain mask on the bold grid
        str14b='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_brain_mask -input brain_mask+tlrc';
        eval(str14b);
        
        %%
        % smoothing is now done in postprocessBold.m
        %         % spatially smooth bolds
        %         str15='!3dBlurInMask -prefix smnudsbold_e1_tlrc_al -input nudsbold_e1_tlrc_al+tlrc -FWHM 5 -mask nudsbold_e1_tlrc_al+tlrc';
        %         eval(str15);
        %
        %         str16='!3dBlurInMask -prefix smnudsbold_e2_tlrc_al -input nudsbold_e2_tlrc_al+tlrc -FWHM 5 -mask nudsbold_e1_tlrc_al+tlrc';
        %         eval(str16);
        %
        %         str17='!3dBlurInMask -prefix smnudsbold_e3_tlrc_al -input nudsbold_e3_tlrc_al+tlrc -FWHM 5 -mask nudsbold_e1_tlrc_al+tlrc';
        %         eval(str17);
        
        cd(origPath);
        
        clear ulf urf urb lrf % needed to determine whether next subject has markers
        clear ulfc urfc urbc lrfc
        
    else % only processing the ROI
        
        if hasMarkers
            cd(pathToData);
            % Afni can't handle overwriting files so...
            delete roi+tlrc.BRIK;
            delete roi+tlrc.HEAD;
            delete controlRoi+tlrc.BRIK;
            delete controlRoi+tlrc.HEAD;
            delete resampled_roi+tlrc.BRIK;
            delete resampled_roi+tlrc.HEAD;
            delete resampled_control_roi+tlrc.BRIK;
            delete resampled_control_roi+tlrc.HEAD;
            
            % align the roi mask to the talairach atlas
            str11='!@auto_tlrc -apar anat+tlrc -input roi+orig';
            eval(str11);
            str11a='!gunzip roi+tlrc.BRIK.gz';
            eval(str11a);
            
            str11b='!@auto_tlrc -apar anat+tlrc -input controlRoi+orig';
            eval(str11b);
            str11c='!gunzip controlRoi+tlrc.BRIK.gz';
            eval(str11c);
            
            % resample the ROI on the bold grid
            %           str13='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_roi -input roi+tlrc';
            str13='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_roi -input roi+tlrc';
            eval(str13);
            
            %             str13b='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_control_roi -input controlRoi+tlrc';
            str13b='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_control_roi -input controlRoi+tlrc';
            eval(str13b);
            cd(origPath);
        else
            % this is done in makeAverageControlROI
            %             %TODO: move this to makeAverageROI
            %             % no markers yet only processing ROI
            %             [err, muRoiMask, InfoMuRoiMask, ~] = BrikLoad (fullfile(pathToData,'muroi+tlrc.BRIK')); % we only do this to get the 'Info' field
            %             % reflect in sagittal direction
            %             muControlRoiMask=muRoiMask;
            %             muControlRoiMask(end:-1:1,:,:)=muRoiMask;
            %             [err,ErrMessage,Info]=WriteBrikWrap(pathToData,muControlRoiMask,InfoMuRoiMask,'muControlRoi','tlrc');
            
        end
        
    end
end