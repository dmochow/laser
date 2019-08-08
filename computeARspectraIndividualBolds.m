% 10.23.18 compute AR spectra and percent change on individual bolds
% NB: significant voxels have been taken from mean-bold analysis
clear all; close all; clc

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

% global parameters
nEchos=3;
ALPHA=0.05;
nTRs=645;
NLAGS=5; %p=NLAGS-1
LASERONSETTR=215;
LASEROFFSETTR=430;

%%
tic
arIndCoeffs=cell(nSubjects,1); %zeros(NLAGS-1,nVoxels,2,nEchos); % no intercept
for s=1:nSubjects
    
    pathToData=['../data/' subjStrs{s} '/NII'];
    brainMaskFilename='resampled_brain_mask+tlrc';
    [~, brainMask, ~, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
    brainMask=logical(brainMask);
    nVoxels=sum(brainMask(:));
    

    pvals=cell(nSubjects,1); %zeros(nVoxels,2);
    tstats=cell(nSubjects,1); %zeros(nVoxels,2);
    
    for e=1:nEchos
        [s,e]
        toc
        dataFilename=['s8oBoldEcho' num2str(e) '+tlrc'];
        [~, bold, info, ~] = BrikLoad (fullfile(pathToData,dataFilename));
  
        boldTs=vol2ts(bold,brainMask);

        for v=1:nVoxels
            %v
            y=boldTs(:,v);
            X=tplitz(y,NLAGS);
            X=X(:,2:end-1); % don't include constant term or present sample
  
            % pre versus stim
            [h1,p1,tstat1] = chowtest(X(1:LASEROFFSETTR,:),y(1:LASEROFFSETTR),LASERONSETTR);
            % pre versus post
            [h2,p2,tstat2] = chowtest(cat(1,X(1:LASERONSETTR,:),X(LASEROFFSETTR+1:end,:)),cat(1,y(1:LASERONSETTR),y(LASEROFFSETTR+1:end)),LASERONSETTR);

            X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
            X2=X(LASERONSETTR+1:LASEROFFSETTR,:); y2=y(LASERONSETTR+1:LASEROFFSETTR);
            X3=X(LASEROFFSETTR+1:end,:); y3=y(LASEROFFSETTR+1:end);
            
            arIndCoeffs{s}(:,v,1,e)=X1\y1;
            arIndCoeffs{s}(:,v,2,e)=X2\y2;
            arIndCoeffs{s}(:,v,3,e)=X3\y3;
            
            pvals{s}(v,1,e)=p1; pvals{s}(v,2,e)=p2;
            tstats{s}(v,1,e)=tstat1;  tstats{s}(v,2,e)=tstat2;
            
        end
        
    end
    
    save('../data/precomputed/arIndividualData.mat','pvals','tstats','arIndCoeffs');
    arIndCoeffs
end




