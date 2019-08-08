clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

% point to the data
path='../data/S10/NII/';
maskFilename='resampled_roi+tlrc';
anatFilename='resampled_anat+tlrc';
epiFilename_1='bold_e1_tlrc_al+tlrc';
epiFilename_2='bold_e2_tlrc_al+tlrc';
epiFilename_3='bold_e3_tlrc_al+tlrc';
TEs=[12.8 34.13 55.46];

% read in the data
[err, mask, Info, ErrMessage] = BrikLoad (fullfile(path,maskFilename));
[err, anat, Info, ErrMessage] = BrikLoad (fullfile(path,anatFilename));
[err, bold_1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
[err, bold_2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
[err, bold_3, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_3));

% compute the mask
finalMask=mask>0 & anat>0; % in the laser beam AND in the brain

% grab the bold in the ROI
tmp=permute(bold_1,[4 1 2 3]);
bold(1,:,:)=tmp(:,:);

tmp=permute(bold_2,[4 1 2 3]);
bold(2,:,:)=tmp(:,:);

tmp=permute(bold_3,[4 1 2 3]);
bold(3,:,:)=tmp(:,:);

%%
bold_roi=bold(:,:,finalMask);

%%

% perform SVD on average of echoes
% TODO: add weighted averaging following Posse et al
% epi1=squeeze(bold_roi(1,:,:));
% epi2=squeeze(bold_roi(2,:,:));
% epi3=squeeze(bold_roi(3,:,:));
muBold=squeeze(mean(bold_roi,1));
[U,S,V]=svd(muBold.',0);
% U are the spatial filters (voxels x components)
% V are associated time courses (time x components)

% normalize time courses
nFrames=size(muBold,1); % number of TRs
nVoxels=size(muBold,2); % number of voxels
nComps=size(U,2); % number of components
tc=V-repmat(mean(V),nFrames,1); % demeaned time courses
tc=tc./repmat(std(V),nFrames,1); % variance normalized time courses
nEchos=3;
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
            psc(c,v,e)=thisBeta(2)/thisBeta(1);
            %psc(c,v,e)=thisBeta(2)/mean(Y);
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
        F1(cIndx,v)=sum((Y-X*fits).^2)./sum(Y.^2);
        
        %
        X=-TEs.';
        fits=(X'*X)\(X'*Y);
        F2(cIndx,v)=sum((Y-X*fits).^2)./sum(Y.^2);
        
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
N=nComps;
sigmas=diag(S);
sigmas=sigmas(sortind);
muBoldClean=U(:,sortind(1:N))*diag(sigmas(1:N))*V(:,sortind(1:N))';

muPre=sum(sum(muBoldClean(:,1:213)));
muStim=sum(sum(muBoldClean(:,214:426)));
muPost=sum(sum(muBoldClean(:,427:645)));

pwrPre=sumsqr(muBoldClean(:,1:213));
pwrStim=sumsqr(muBoldClean(:,214:426));
pwrPost=sumsqr(muBoldClean(:,427:645));

% pwrPre=sum(sum(muBold(1:213,:)));
% pwrStim=sum(sum(muBold(214:426,:)));
% pwrPost=sum(sum(muBold(427:645,:)));

figure;
subplot(131);
plot(mean(muBoldClean,1))
subplot(132);
bar([muPre muStim muPost]);
subplot(133);
bar([pwrPre pwrStim pwrPost]);


% meantc=mean(tc(:,d(1:200)),2);
% figure; plot(meantc); 
% 
% sum(meantc(1:213))
% sum(meantc(213:426))
% sum(meantc(426:end))

% sum(muBold(1:213,4))
% sum(muBold(213:426,4))





    


%% assume biopac in memory as data
% [~,maxindsample]=max(diff(data(:,1)));
% onsetTimeSec=maxindsample/1000;
% onsetTimeTR=onsetTimeSec/2.8;

