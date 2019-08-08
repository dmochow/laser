%% 05.05.19 analysis on MEICA preprocessed data
clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
%subjStrs={'S23','S24','S25','S26','S27','S28','S29','S30','S31','S32'};
subjStrs={'S23'};
nSubjects=numel(subjStrs);
alpha=0.05;
nTRsMax=685; % nTRs in the longest recording
TR=2.8;
anatRoiMaskFilename='roi_r21_z39+tlrc.BRIK';
anatBrainMaskFilename='brain_mask+tlrc.BRIK';
precomputedFilename='../data/precomputed/ROIGLMToggleStats.mat';
roiMaskFilename='follow_ROI_LASER+tlrc.BRIK';
segFilename='resampled_Classes+tlrc.BRIK';
gmIndex=2;
nTRsRemoved=3; % MEICA

for s=1:nSubjects
    s
    pathToData=['../data/output/' subjStrs{s} '/' subjStrs{s} '.results'];
    pathToBiopac=['../data/scanner_data/' subjStrs{s}];
    boldFilename=['all_runs.' subjStrs{s} '+tlrc.BRIK'];
    
    [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,boldFilename));
      
     % illumination ROI
    [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
    roiMask=logical(roiMask);
    
    % grey matter mask
    [~, seg, ~, ~] = BrikLoad (fullfile(pathToData,segFilename));
    greyMask=logical(seg==gmIndex);
    
    finalMask=roiMask&greyMask;

    
    %%
    biopacFilename=dir(fullfile(pathToBiopac,'subj*.mat'));
    load(fullfile(pathToBiopac,biopacFilename.name),'data');
    thresh=0.75; % anything over this voltage is classified "laser on"
    laserTc=data;
    biopacTime=(0:size(data,1)-1)*0.001; % biopac time in ms but check by loading in all the variables
    boldTime=(0:644)*TR; % note that the biopac time is likely shorter than bold
    laserTcDown=interp1(biopacTime,laserTc(:,1),boldTime);
    regressor=laserTcDown>thresh;
    
    %%
    nTRs=size(bold,4)+nTRsRemoved; % 
    regressor=cat(2,regressor,zeros(1,nTRs-numel(regressor)));
    regressor=regressor(nTRsRemoved+1:end);
    X=regressor(:);
 
%     if numel(X) < nTRs
%         X=cat(1,X,zeros(nTRs-nn,1));
%     end
    Xnull=ones(size(X));
    allRegressors{1,s}=X;
    
    %%
    ts=vol2ts(bold,finalMask);
    
%     % fill out ts to max number of TRs
%     nn=size(ts,1);
%     if nn < nTRs
%         ts=cat(1,ts,zeros(nTRs-nn,size(ts,2)));
%     end
    
    allts{s}=ts;
    nVoxels=size(ts,2);
    for v=1:nVoxels
        statsf=myGLM(ts(:,v),X,1); % full model
        statsr=myGLM(ts(:,v),Xnull,0); % reduced model
        allBetas{s}(v)=statsf.B(1);
        ssef=statsf.sse;
        doff=statsf.dof;
        sser=statsr.sse;
        dofr=statsr.dof;
        Fstat=((sser-ssef)/(dofr-doff)) / (ssef/doff);
        pval{s}(v)=fcdf(Fstat,dofr-doff,doff,'upper');
    end
    
    [p_fdr,p_masked]=fdr(pval{s},alpha);
    issig{s}=p_masked;
    
%     % convert issig to afni brik and resample to anatomical space
%     tmp=issig{e,s};
%     issigvol=ts2vol(tmp,finalMask); % issig in bold space
%     [err,ErrMessage,Info]=WriteBrikWrap(pathToData,issigvol,boldInfo,['issigEcho' num2str(e)],'tlrc');
%     
%     % convert betas to afni brik and resample to anatomical space
%     tmp=allBetas{e,s};
%     betavol=ts2vol(tmp,finalMask); % betas in bold space
%     [err,ErrMessage,Info]=WriteBrikWrap(pathToData,betavol,boldInfo,['toggleBetasEcho' num2str(e)],'tlrc');
%     
%     % now resample to anat space
%     origPath=pwd; cd(pathToData);
%     prefixStr=['issigAnatEcho' num2str(e)];
%     inputStr=['issigEcho' num2str(e) '+tlrc'];
%     afniStr=['!3dresample -master anat+tlrc -prefix ' prefixStr ' -input ' inputStr]; eval(afniStr);
%     
%     prefixStr=['toggleBetasAnatEcho' num2str(e)];
%     inputStr=['toggleBetasEcho' num2str(e) '+tlrc'];
%     afniStr=['!3dresample -master anat+tlrc -prefix ' prefixStr ' -input ' inputStr]; eval(afniStr);
%     
%     cd(origPath);
    
end



%%
r=cellfun(@mean,issig);
r2=cellfun(@sum,issig);

%%
%save and exit
save(precomputedFilename,'issig','allBetas','allts','allRegressors','r','r2');

%%
s=1; sigBetas=allBetas{s}(issig{s});

%%
% look at the significant voxels
% duration=round(180/TR)-2;
% figure; hold on
% for e=1:3
%     for s=1
%         sigvox=find(issig{e,s});
%         tcs=allts{e,s}(:,sigvox);
%         %plot(mean(tcs,2)); %plot(regressor);
%         ontimes=find(diff(regressor)==1);
%         epochs = simpleEpoch(tcs',ontimes,duration);
%         muEpoch=mean(epochs,3);
%         muMuEpoch=mean(muEpoch,1);
%         %hs(e,s)=subplot(3,nSubjects,(e-1)*nSubjects+s);
%         hs(e,s)=subplot(3,nSubjects,e);
%         plot(muMuEpoch);
%     end
% end

% %%
% allsigts=cell(3,1);
% for e=1:3
%     for s=1:nSubjects
%         tmp=allts{e,s}(:,issig{e,s});
%         nn=size(tmp,1);
%         if nn < nTRs
%             tmp=cat(1,tmp,zeros(nTRs-nn,size(tmp,2)));
%         end
%         allsigts{e}=cat(2,allsigts{e},tmp);
%     end
% end
%
% %%
% figure;
% for e=1:3
%     subplot(3,1,e);
%     plot(mean(allsigts{e},2));
% end

%%
% % figure
% hf=figure;
% nRows=5;
% % %significant vs subject
% for s=1:nSubjects
%     
%     zcoor=75;
%     pathToData=['../data/' subjStrs{s}];
%     [~, anat, ~, ~] = BrikLoad (fullfile(pathToData,'anat+tlrc'));
%     anat=anat(end:-1:1,:,:);
%     
%     [~, anatRoiMask, ~, ~] = BrikLoad (fullfile(pathToData,anatRoiMaskFilename));
%     anatRoiMask=anatRoiMask(end:-1:1,:,:);
%     anatRoiMask=logical(anatRoiMask);
%     
%     [~, anatBrainMask, ~, ~] = BrikLoad (fullfile(pathToData,anatBrainMaskFilename));
%     anatBrainMask=anatBrainMask(end:-1:1,:,:);
%     
%     [~, issigMask, ~, ~] = BrikLoad (fullfile(pathToData,['issigAnatEcho1+tlrc']));
%     issigMask=issigMask(end:-1:1,:,:);
%     
%     
%     hs(1,s)=subplot(nRows,nSubjects,s);
%     h = colorBrain(anat(:,:,zcoor),anatRoiMask(:,:,zcoor),issigMask(:,:,zcoor)==1,anatBrainMask(:,:,zcoor),1);
%     
%     %        hbar=bar(1:nSubjects,100*r');
%     %        set(get(hs(1),'ylabel'),'String','% significant voxels');
%     %        set(get(hs(1),'xlabel'),'String','Subject index');
%     
% end

% % two time courses
% s=4;
% hs(2)=subplot(5,2,3); hold on
% for e=1:3
%     sigvox=find(issig{e,s});
%     tcs=allts{e,s}(:,sigvox);
%     plot(mean(tcs,2));
%     plot(regressor*0.5,'k');
% end
%
% s=6;
% hs(3)=subplot(5,2,4); hold on
% for e=1:3
%     sigvox=find(issig{e,s});
%     tcs=allts{e,s}(:,sigvox);
%     plot(mean(tcs,2));
%     plot(regressor*0.5,'k');
% end
%
% print -dpng ../figures/tmp


