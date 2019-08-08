clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
subjStrs={'S23','S24','S25','S26','S27','S28','S29','S30','S31','S32'};
nSubjects=numel(subjStrs);
alpha=0.05;
%nTRs=685; % nTRs in the longest recording
TR=2.8;
anatRoiMaskFilename='roi_r21_z39+tlrc.BRIK';
anatBrainMaskFilename='brain_mask+tlrc.BRIK';
savedDataFilename='../data/precomputed/revFig2data.mat';
gmIndex=2;

for e=1:3
    echoStr=num2str(e);
    filename=['smoothedOutputBoldEcho' echoStr '+tlrc.BRIK'];
    
    for s=1:nSubjects
        [e,s]
        pathToData=['../data/myoutput/' subjStrs{s}]
        
        [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,filename));
        nTRs=size(bold,4);
        
        brainMaskFilename='resampled_brain_mask+tlrc';
        [~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
        brainMask=logical(brainMask);
        
        %grey matter mask
        resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
        [~, seg, info, ~] = BrikLoad (resampledSegFilename);
        brainMask=logical(seg==gmIndex);
%         
        roiMaskFilename='resampled_roi_r21_z39+tlrc.BRIK';
        [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
        roiMask=logical(roiMask);
        
        finalMask=roiMask&brainMask;
        
        load(fullfile(pathToData,'regressor'));
        regressor=cat(2,regressor,zeros(1,nTRs-numel(regressor)));
        X=regressor(:);
        if numel(X) < nTRs
            X=cat(1,X,zeros(nTRs-nn,1));
        end
        Xnull=ones(size(X));
        allRegressors{1,s}=X;
        
        %%
        ts=vol2ts(bold,finalMask);
        
        % fill out ts to max number of TRs
%         nn=size(ts,1);
%         if nn < nTRs
%             ts=cat(1,ts,zeros(nTRs-nn,size(ts,2)));
%         end
        
        allts{e,s}=ts;
        nVoxels=size(ts,2);
        for v=1:nVoxels
            statsf=myGLM(ts(:,v),X,1); % full model
            statsr=myGLM(ts(:,v),Xnull,0); % reduced model
            allBetas{e,s}(:,v)=statsf.B;
            ssef=statsf.sse;
            doff=statsf.dof;
            sser=statsr.sse;
            dofr=statsr.dof;
            Fstat=((sser-ssef)/(dofr-doff)) / (ssef/doff);
            pval{s}(v)=fcdf(Fstat,dofr-doff,doff,'upper');
        end
        
        [p_fdr,p_masked]=fdr(pval{s},alpha);
        issig{e,s}=p_masked;
        
%         % convert issig to afni brik and resample to anatomical space
%         tmp=issig{e,s};
%         issigvol=ts2vol(tmp,finalMask); % issig in bold space
%         [err,ErrMessage,Info]=WriteBrikWrap(pathToData,issigvol,boldInfo,['issigEcho' num2str(e)],'tlrc');
%         
%         % convert betas to afni brik and resample to anatomical space
%         tmp=allBetas{e,s};
%         betavol=ts2vol(tmp,finalMask); % betas in bold space
%         [err,ErrMessage,Info]=WriteBrikWrap(pathToData,betavol,boldInfo,['toggleBetasEcho' num2str(e)],'tlrc');
% 
%         % now resample to anat space
%         origPath=pwd; cd(pathToData);
%         prefixStr=['issigAnatEcho' num2str(e)];
%         inputStr=['issigEcho' num2str(e) '+tlrc'];
%         afniStr=['!3dresample -master anat+tlrc -prefix ' prefixStr ' -input ' inputStr]; eval(afniStr);
%         
%         prefixStr=['toggleBetasAnatEcho' num2str(e)];
%         inputStr=['toggleBetasEcho' num2str(e) '+tlrc'];
%         afniStr=['!3dresample -master anat+tlrc -prefix ' prefixStr ' -input ' inputStr]; eval(afniStr);
%         
%         cd(origPath);

    end
    
end

%%
r=cellfun(@mean,issig);
r2=cellfun(@sum,issig);

%%
% save and exit
save(savedDataFilename,'issig','allBetas','allts','allRegressors','r','r2');


