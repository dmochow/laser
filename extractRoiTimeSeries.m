clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

% point to the data
path='../data/S02/NII/';
maskFilename='resampled_roi+tlrc';
altMaskFilename='muroi+tlrc';
anatFilename='resampled_anat+tlrc';
epiFilename_1='smdsbold_e1_tlrc_al+tlrc';
epiFilename_2='smdsbold_e2_tlrc_al+tlrc';
epiFilename_3='smdsbold_e3_tlrc_al+tlrc';
TEs=[12.8 34.13 55.46];
nEchos=numel(TEs);
% TODO: import biopac information here
% ask for path, then compute laser onset/offset time

% read in the data
[err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,maskFilename));
if err
    [err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,altMaskFilename));
end
[err, anat, Info, ErrMessage] = BrikLoad (fullfile(path,anatFilename));
[err, bold_1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
[err, bold_2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
[err, bold_3, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_3));

% compute the mask
finalMask=mask>0 & anat>0; % in the laser beam AND in the brain
numel(find(finalMask))
% grab the bold in the ROI
tmp=permute(bold_1,[4 1 2 3]);
bold(1,:,:)=tmp(:,:);

tmp=permute(bold_2,[4 1 2 3]);
bold(2,:,:)=tmp(:,:);

tmp=permute(bold_3,[4 1 2 3]);
bold(3,:,:)=tmp(:,:);

%%
bold_roi=bold(:,:,finalMask);
% echo x TR x voxel

%%
% perform SVD on average of echoes
tmp=permute(bold_roi,[2 3 1]);
catBold=tmp(:,:);

% preprocess for PCA
catBold_dem=catBold-repmat(mean(catBold,2),1,size(catBold,2));
catBold_varnorm=catBold./repmat(std(catBold,[],2),1,size(catBold,2));
% PCA...
[U,S,V]=svd(catBold_varnorm,0);
nPCs2Keep=75; %TODO: implement an automatic selection rule
diagS=diag(S);
catBoldTrunc=U(:,1:nPCs2Keep)*diag(diagS(1:nPCs2Keep))*V(:,1:nPCs2Keep).';

% ...then ICA
[icasig, A, W]=fastica(catBoldTrunc.');  % with transpose: voxels x TRs
% with the transpose, this treats each voxel as a sensor with a time course
tc=icasig.'; % each column is time course

% normalize time courses
nFrames=size(catBoldTrunc,1); % number of TRs
nVoxels=size(catBoldTrunc,2)/nEchos; % number of voxels
nComps=size(tc,2); % number of components
tc=tc-repmat(mean(tc),nFrames,1); % demeaned time courses
tc=tc./repmat(std(tc),nFrames,1); % variance normalized time courses

%%
% compute betas
% regress bold onto time course of each component at each voxel and echo
for c=1:nComps
    c
    for v=1:nVoxels
        for e=1:nEchos
            Y=squeeze(bold_roi(e,:,v)).';
            X=[ ones(nFrames,1) tc(:,c)  ];
            thisBeta=(X'*X)\(X.'*Y);
            beta1(c,v,e)=thisBeta(2);
            %psc(c,v,e)=thisBeta(2)/thisBeta(1);
            psc(c,v,e)=thisBeta(2)/mean(Y);
        end
    end
end

% components x voxels
%%
F1=zeros(nComps,nVoxels);
F2=zeros(nComps,nVoxels);
for cIndx=1:nComps
    for v=1:nVoxels
        Y=squeeze(psc(cIndx,v,:));
        
        %
        X=ones(3,1);
        fits=(X'*X)\(X'*Y);
        %F1(cIndx,v)=sum((Y-X*fits).^2)./sum(Y.^2);
        F1(cIndx,v)=(sum(Y.^2)-sum((Y-X*fits).^2))/sum(Y.^2);
        
        %
        X=TEs.';
        fits=(X'*X)\(X'*Y);
        F2(cIndx,v)=(sum(Y.^2)-sum((Y-X*fits).^2))/sum(Y.^2);
        
    end
end

alpha=sum(beta1.^2,3);
kappa=sum(alpha.*F2,2)./sum(alpha,2);
rho=sum(alpha.*F1,2)./sum(alpha,2);

figure;
subplot(211);
plot(sort(kappa,'descend'),'*');
subplot(212);
plot(sort(rho,'descend'),'*');


%%
[~,sortind]=sort(kappa,'descend');
%N=nPCs2Keep;
N=2;
%inds2keep=sortind(1:50);
inds2keep=sortind(1:N);
%inds2keep=sortind(end-N+1:1:end);
muBoldClean=(A(:,inds2keep)*icasig(inds2keep,:)).';

muPre=sum(sum(muBoldClean(1:213,:)));
muStim=sum(sum(muBoldClean(214:426,:)));
muPost=sum(sum(muBoldClean(427:645,:)));

pwrPre=sumsqr(muBoldClean(1:213,:));
pwrStim=sumsqr(muBoldClean(214:426,:));
pwrPost=sumsqr(muBoldClean(427:645,:));

%
x=mean(muBoldClean,2);
%y = detrend(x,'linear',[1 numel(x)]);

%
TR=2.8;
figure;
subplot(211);
plot((1:645)*TR,mean(muBoldClean,2))
%plot((1:645)*TR,y)
subplot(223);
bar([muPre muStim muPost]);
subplot(224);
bar([pwrPre pwrStim pwrPost]);

%print -dpng ../figures/jigga_very_early
%% assume biopac in memory as data
% [~,maxindsample]=max(diff(data(:,1)));
% onsetTimeSec=maxindsample/1000;
% onsetTimeTR=onsetTimeSec/2.8;


