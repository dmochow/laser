function cshifts=computeCshifts(option)
% compute the spatial shifts to align everybody's brain in the talairach
% space
% note that this takes everyone away, just slightly, from talairach

pathToData='../data/SAVG/NII';
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17',...
    'S18','S19','S20','S21','S22'}; % all subjects with markers
nSubjects=numel(subjStrs);

switch option
    case 1 % use the laser origin (doesn't account for orientation)
        
        %% need laserOrigin, urf,... so...
        % voxel space to mri space
        % mri space to talairach mri space
        % talairach mri space to voxel space
        XYZOUT=zeros(3,nSubjects);
        for s=1:nSubjects
            subjStr=subjStrs{s};
            pathToData=['../data/' subjStr '/NII/'];
            
            % get laser origin
            try
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
                xyzOut = orig2tal(Xin,Xat,InfoOrig,InfoTlrc);
            catch
                xyzOut=NaN(3,1);
            end
            XYZOUT(:,s)=xyzOut;
        end
        muXYZ=nanmean(XYZOUT,2);
        cshifts_xyz=XYZOUT-repmat(muXYZ,1,nSubjects);
        sizeRatio=2.5/0.9;
        cshifts_xyz_resampled=cshifts_xyz/sizeRatio;
        cshifts=round(cshifts_xyz_resampled);
        cshifts(isnan(cshifts))=0; % replance missing data with "no shift"
        
    case 2 % using centroid of "ROI"
        nX=64; nY=76; nZ=60; % dimensions of BOLD
        radMask=19.2; depthMask=26.5; % original values
        roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
        muRoiPrefix=['mu_' roiPrefix];
        centroid=zeros(3,nSubjects);
        for s=1:nSubjects
            subjIndx=subjStrs{s};
            basePath=['../data/' subjIndx '/NII/'];
            brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');
            roiMaskFilename=fullfile(basePath,['resampled_' roiPrefix '+tlrc.BRIK']);
            muRoiMaskFilename=fullfile(basePath,['resampled_' muRoiPrefix '+tlrc.BRIK']);
            [~, roiMask, ~, ~] = BrikLoad (roiMaskFilename);
            if isempty(roiMask)
                [~, roiMask, ~, ~] = BrikLoad (muRoiMaskFilename);
            end
            roiMask=logical(roiMask);
            [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
            brainMask=logical(brainMask);
            finalMask=roiMask>0 & brainMask>0;
            [X,Y,Z]=ndgrid(1:nX,1:nY,1:nZ);
            xmask=X(finalMask); ymask=Y(finalMask); zmask=Z(finalMask);
            centroid(:,s)=[mean(xmask);mean(ymask);mean(zmask)];
        end
        muCentroid=mean(centroid,2); % try to align everyone's bold to for a coherent spatial average
        cshifts=round(centroid-repmat(muCentroid,1,nSubjects));
end

