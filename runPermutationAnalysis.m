clear all; close all; clc

% permutation analysis to detect significant differences in group-mean
% waveform
dataFilename='../data/precomputed/subjectLevelGlmResults 20-Mar-2018';
load(dataFilename);
nPerms=300;
nSubjects=numel(allBoldsRoi);
nEchos=3;
%%
allBoldsMat=cat(3,allBoldsRoi{:});
muBolds=mean(allBoldsMat,3);
for e=1:nEchos
    mockBold=[];
    for s=1:nSubjects
        thisBold=squeeze(allBoldsRoi{s}(e,:,:));
        [nFrames,nVoxels]=size(thisBold);
        thisMockBold=zeros(nFrames,nVoxels,nPerms);
        for p=1:nPerms
            [e,s,p/nPerms*100]
            thisMockBold(:,:,p)=surrogateResponseGenerator(thisBold);
        end
        mockBold=cat(2,mockBold,thisMockBold);
    end
    
    muMockBold=squeeze(mean(mockBold,2));
    % at each time point, evaluate significant deviation of true record
    pvals(:,e)=mean(muMockBold>repmat(muBolds(e,:).',[1 nPerms]),2);
    
%     
%     figure; hold on
%     plot(muMockBolds);
%     plot(muBolds,'k','LineWidth',4);
    
end

%%
alpha=0.05;
[p_fdr, p_masked] = fdr( pvals, alpha);

%%
TR=2.8;
time=(0:nFrames-1)*TR;
figure;
for e=1:nEchos
    subplot(nEchos,1,e); hold on;
    plot(time,squeeze(mean(allBoldsMat(e,:,:),3)));
    stem(time(find(p_masked(:,e))),zeros(sum(p_masked(:,e),1)),'*');
end

%%
for p=1:nPerms
    mockMus(:,:,p)=surrogateResponseGenerator(muBolds.');
end

pvals=mean(mockMus>repmat(muBolds.',[1 1 nPerms]),3);

figure; 
plot(squeeze(mockMus(:,1,:)));
hold on
plot(muBolds(1,:),'k','LineWidth',5);

