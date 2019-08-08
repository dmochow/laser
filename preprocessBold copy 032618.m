clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

% we need the following files (and only those files) present in the pathToData folder:
% anat.nii
% bold_e1.nii
% bold_e2.nii
% bold_e3.nii

subjStrs={'S04'};
bigmove=0; % parameter for AfNI alignment (have not seen it improve anything)
%subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'};
nSubjects=numel(subjStrs);

% important parameters
RADMASK=19.2; % lateral laser reach
DEPTHMASK=26.5; % depth laser reach

for s=1:nSubjects
    subjStr=subjStrs{s};
    
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
    end
    pathToData=['../data/' subjStr '/NII/'];
    
    % compute brain mask in afni
    origPath=pwd;
    cd(pathToData);
    if exist('brain_mask+orig.BRIK','file')
        delete brain_mask+orig.BRIK;
        delete brain_mask+orig.HEAD;
    end
    str1='!3dcopy anat.nii anat';
    str2='!3dSkullStrip -input anat+orig';
    str3='!3dAutomask -prefix brain_mask skull_strip_out+orig';
    eval(str1);
    eval(str2);
    eval(str3);
    cd(origPath);
    
    % processing begins here
    pathToData=['../data/' subjStr '/NII/'];
    niiFilename=['../data/' subjStr '/NII/anat.nii'];
    anatFilename=['../data/' subjStr '/NII/anat+orig'];
    brainMaskFilename=['../data/' subjStr '/NII/brain_mask+orig'];
    outNiifilename=['../data/' subjStr '/NII/roiMask.nii'];
    outControlNiifilename=['../data/' subjStr '/NII/controlRoiMask.nii'];
    
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
        controlLaserOrigin = getLaserOrigin(ulfc,urfc,urbc);
        
        % fine-tune the position of real laser (project it onto the scalp)
        laserOrigin = fineTuneOrigin(anatFilename,laserOrigin,ulf,urf,lrf);
        save(fullfile(pathToData,'laserOrigin.mat'),'laserOrigin');
        
        % fine-tune the position of the control laser
        controlLaserOrigin = fineTuneOrigin(anatFilename,controlLaserOrigin,ulfc,urfc,lrfc);
        if isempty(controlLaserOrigin) % this happens sometimes
            % solution is to reflect the fine-tuned true laser
            controlLaserOrigin=laserOrigin; % AFTER fine-tuning
            controlLaserOrigin(1)=2*midsag-laserOrigin(1);
        end
        save(fullfile(pathToData,'controlLaserOrigin.mat'),'controlLaserOrigin');
        
        % project the real laser into the mri
        [err, mask, Info, ErrMessage] = BrikLoad (brainMaskFilename);
        mridim=size(mask);
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
        roiMask=zeros(size(mask));
        roiMask( (xr.^2+zr.^2)<RADMASK.^2 & yr>0 & yr<DEPTHMASK )=1;
        
        %% write out ROI mask
        InfoWrite=Info;
        opt.Scale=0;
        opt.Prefix='roi';
        opt.View='orig';
        opt.Verbose=1;
        opt.AppendHistory=1;
        opt.NoCheck=0;
        opt.Overwrite=1;
        origPath=pwd;
        cd(pathToData);
        if exist('roi+orig.BRIK','file')
            delete roi+orig.BRIK;
            delete roi+orig.HEAD;
        end
        [err, ErrMessage, Info] = WriteBrik (roiMask, InfoWrite,opt);
        cd(origPath);
        
    end
    
    % project the control laser onto the mri
    cLOx=controlLaserOrigin(1)*ones(size(x));
    cLOy=controlLaserOrigin(2)*ones(size(y));
    cLOz=controlLaserOrigin(3)*ones(size(z));
    
    %%
    % translation
    T=eye(4);
    T(1,4)=-controlLaserOrigin(1);
    T(2,4)=-controlLaserOrigin(2);
    T(3,4)=-controlLaserOrigin(3);
    
    xt=T(1,1)*x+T(1,2)*y+T(1,3)*z+T(1,4)*1;
    yt=T(2,1)*x+T(2,2)*y+T(2,3)*z+T(2,4)*1;
    zt=T(3,1)*x+T(3,2)*y+T(3,3)*z+T(3,4)*1;
    
    %%
    % rotation
    u=(urfc-ulfc)/norm(urfc-ulfc);  % LR
    v=(urfc-urbc)/norm(urfc-urbc);  % PA
    w=(urfc-lrfc)/norm(urfc-lrfc);  % IS
    R=eye(4);
    R(1:3,1)=u;
    R(1:3,2)=v;
    R(1:3,3)=w;
    R=inv(R);
    
    xr=R(1,1)*xt+R(1,2)*yt+R(1,3)*zt+R(1,4)*1;
    yr=R(2,1)*xt+R(2,2)*yt+R(2,3)*zt+R(2,4)*1;
    zr=R(3,1)*xt+R(3,2)*yt+R(3,3)*zt+R(3,4)*1;
    
    %%
    controlRoiMask=zeros(size(mask));
    controlRoiMask( (xr.^2+zr.^2)<RADMASK.^2 & yr>0 & yr<DEPTHMASK )=1;
    
    %% write out control ROI mask
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
    [err, ErrMessage, Info] = WriteBrik (controlRoiMask, InfoWrite,opt);
    cd(origPath);
    
    
    
    
    %%
    % entering what used to be 'runAfniPreProcessing.m'
    cd(pathToData);
    str4='!3dcopy bold_e1.nii bold_e1';
    eval(str4);
    
    str5='!3dcopy bold_e2.nii bold_e2';
    eval(str5);
    
    str6='!3dcopy bold_e3.nii bold_e3';
    eval(str6);
    
    str9='!3ddespike -NEW -localedit -prefix dsbold_e1 bold_e1+orig';
    eval(str9);
    
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
    
    if hasMarkers
        % align the roi mask to the talairach atlas
        str11='!@auto_tlrc -apar anat+tlrc -input roi+orig';
        eval(str11);
        
        str11b='!@auto_tlrc -apar anat+tlrc -input controlRoi+orig';
        eval(str11b);
    end
    
    % now align the epis to the anatomy, and also telling it that we have an
    % atlased anatomy
    
    if ~bigmove
        str12='!align_epi_anat.py -anat anat+orig -epi dsbold_e2+orig -epi_base 5 -child_epi dsbold_e1+orig dsbold_e3+orig -epi2anat -suffix _al -tlrc_apar anat+tlrc';
        eval(str12);
    else
        str12='!align_epi_anat.py -big_move -anat anat+orig -epi dsbold_e2+orig -epi_base 5 -child_epi dsbold_e1+orig dsbold_e3+orig -epi2anat -suffix _al -tlrc_apar anat+tlrc';
        eval(str12);
    end
    
    str12a='!3dTproject -input dsbold_e1_tlrc_al+tlrc -prefix nudsbold_e1_tlrc_al -ort dsbold_e1_vr_motion.1D';
    eval(str12a);
    
    str12b='!3dTproject -input dsbold_e2_tlrc_al+tlrc -prefix nudsbold_e2_tlrc_al -ort dsbold_e2_vr_motion.1D';
    eval(str12b);
    
    str12c='!3dTproject -input dsbold_e3_tlrc_al+tlrc -prefix nudsbold_e3_tlrc_al -ort dsbold_e3_vr_motion.1D';
    eval(str12c);
    
    
    if hasMarkers
        % resample the ROI on the bold grid
        str13='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_roi -input roi+tlrc';
        eval(str13);
        
        str13b='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_control_roi -input controlRoi+tlrc';
        eval(str13b);
    end
    
    % resample the atlased anatomy on the bold grid?
    str14='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_anat -input anat+tlrc';
    eval(str14);
    
    % resample the brain mask on the bold grid
    str14b='!3dresample -master nudsbold_e2_tlrc_al+tlrc -prefix resampled_brain_mask -input brain_mask+tlrc';
    eval(str14b);
    
    %%
    % spatially smooth bolds
    str15='!3dBlurInMask -prefix smnudsbold_e1_tlrc_al -input nudsbold_e1_tlrc_al+tlrc -FWHM 5 -mask nudsbold_e1_tlrc_al+tlrc';
    eval(str15);
    
    str16='!3dBlurInMask -prefix smnudsbold_e2_tlrc_al -input nudsbold_e2_tlrc_al+tlrc -FWHM 5 -mask nudsbold_e1_tlrc_al+tlrc';
    eval(str16);
    
    str17='!3dBlurInMask -prefix smnudsbold_e3_tlrc_al -input nudsbold_e3_tlrc_al+tlrc -FWHM 5 -mask nudsbold_e1_tlrc_al+tlrc';
    eval(str17);
    
    cd(origPath);
    
    clear ulf urf urb lrf % needed to determine whether next subject has markers
    clear ulfc urfc urbc lrfc
    
end