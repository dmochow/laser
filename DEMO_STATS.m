clear all; close all; clc
pathToData='.';
bootFilenames={'boostrappedGroupStatsEcho1.mat','boostrappedGroupStatsEcho2.mat','boostrappedGroupStatsEcho3.mat'};

% load echo 1 statistics
load(fullfile('../data/precomputed',bootFilenames{1}),'sigTs','nsigTs','sigMuSpectra','nsigMuSpectra');
TR=2.8;
fs=1/TR;
freqs=[0:size(sigMuSpectra,1)-1]/size(sigMuSpectra,1)*fs; % frequency in Hz
muSpectra=nanmean(sigMuSpectra,3); % means
stdSpectra=nanstd(sigMuSpectra,[],3); % standard deviations (not plotting these below, just computing them here)

% draw
figure; 
hstem=stem(freqs,[sigMuSpectra(:,1) sigMuSpectra(:,2) sigMuSpectra(:,3)],'filled');
set(hstem(3),'MarkerEdgeColor',[0.7 0.7 0.7]);
set(hstem(3),'MarkerFaceColor',[0.7 0.7 0.7]);
xlim([0 fs/2])


%%
% revision because NKI has no matlab
for b=1:3
    load(fullfile('../data/precomputed',bootFilenames{b}),'sigMuSpectra','nsigMuSpectra');
    vals1=permute(sigMuSpectra,[1 3 2]);
    csvwrite(['../data/precomputed/statsBeforeLaserEcho' num2str(b) '.csv'],vals1(:,:,1));
    csvwrite(['../data/precomputed/statsDuringLaserEcho' num2str(b) '.csv'],vals1(:,:,2));
    csvwrite(['../data/precomputed/statsAfterLaserEcho' num2str(b) '.csv'],vals1(:,:,3));
    
    vals2=permute(nsigMuSpectra,[1 3 2]);
    csvwrite(['../data/precomputed/nullStatsBeforeLaserEcho' num2str(b) '.csv'],vals2(:,:,1));
    csvwrite(['../data/precomputed/nullStatsDuringLaserEcho' num2str(b) '.csv'],vals2(:,:,2));
    csvwrite(['../data/precomputed/nullStatsAfterLaserEcho' num2str(b) '.csv'],vals2(:,:,3));
end





