clear all; close all; clc

load ../data/precomputed/groupRoiData9subjects allBolds allOnsets allMaskInds

%%
catBolds=cat(3,allBolds{:});
muCatBold=squeeze(mean(catBolds,1));
%%
% PCA...
[U,S,V]=svd(muCatBold,0);
nPCs2Keep=150; %TODO: implement an automatic selection rule
diagS=diag(S);
muBoldTrunc=U(:,1:nPCs2Keep)*diag(diagS(1:nPCs2Keep))*V(:,1:nPCs2Keep).';

% ...then ICA
[icasig, A, W]=fastica(muBoldTrunc.');  % with transpose: voxels x TRs
% with the transpose, this treats each voxel as a sensor with a time course
tc=icasig.'; % each column is time course