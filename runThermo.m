clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%% CONSTANTS
a=-0.01e-6; % "a" Yuan et al. 2012
gamma=2.675e8; % gyromagnetic ratio of hydrogen in water
Bo=3; % 3 Tesla
TE=17e-3; % 17 ms?
delrad=8*pi/(2^16-1); %
gmIndex=2;

%%
%10-minute scan
% nTRs=150;
% TR=20.5*60/nTRs; % HARD-CODED FOR NOW
% time=(0:nTRs-1)*TR;
% onsetTR=21; %
% offsetTR=21+round(10*60/TR);
% subjStrs={'tS01','tS01b','tS01c'};
% %subjStrs={'tS01b','tS01c'};
% %titlestr={'AJ 07.31.19'};
% titles={'AJ July 18 laser on'; 'AJ July 31 laser off'; 'AJ July 31 laser on' };

% 5-minute scan
nTRs=50;
TR=7*60/nTRs; % HARD-CODED FOR NOW
time=(0:nTRs-1)*TR;
onsetTR=6; % 
offsetTR=43;
subjStrs={'tS02b','tS02a'};
%titlestr={'JD 07.31.19'};
titles={'JD July 31 laser off'; 'JD July 31 laser on'};
%%

slope=8*pi/2^16;
intercept=-4*pi;
P=4; Q=4; zero=1; override=1; % params for bad voxels
nPCs=3;
%
nSubj=numel(subjStrs);

for s=1:nSubj
    
    subjStr=subjStrs{s};
    pathToData=['../data/thermo/output/' subjStr '/'];
    
    inputPhFilename=fullfile(pathToData,['ph+orig.BRIK']);
    [~, inputPh, iInfo, ~] = BrikLoad (inputPhFilename);
    
    tempChange=rem(inputPh,pi); % (-4*pi,4*pi) --> (-pi,pi)
    
    tempChange=tempChange/(a*gamma*Bo*TE); % PRF formula
       
    roiMaskFilename=fullfile(pathToData,['resampled_roi_r21_z39+orig.BRIK']);
    [~, roiMask, ~ , ~] = BrikLoad (roiMaskFilename);
    
    % grey matter mask
    resampledSegFilename=fullfile(pathToData,'resampled_Classes+orig.BRIK');
    [~, seg, ~, ~] = BrikLoad (resampledSegFilename);
    greyMask=logical(seg==gmIndex);
    
    finalMask=roiMask&greyMask;
    roiTempChange=vol2ts(tempChange,finalMask);
    
    % bring in drift file
    fid=fopen(fullfile(pathToData,'drift.txt'));
    if fid>0
        formatSpec=['rep %u phase = %f \n'];
        fgets(fid); fgets(fid);
        D = fscanf(fid,formatSpec);
        fclose(fid);
        drift=D(2:2:end);   
        % regress out drift
        roiTempChange=regressOut(roiTempChange.',drift',1).';
    end
    
    %%
    regMask=~roiMask&greyMask;
    regTemp=vol2ts(tempChange,regMask);
    [U,S,V]=svd(regTemp,0);
    roiTempClean=regressOut(roiTempChange.',U(:,1:nPCs).',1).';
%     roiTempClean=roiTempChange;
    
    % remove bad voxels
    %roiTempChange=nanBadChannels(roiTempChange.',P,Q,zero,override).';
    roiTempClean=nanBadChannels(roiTempClean.',P,Q,zero,override).';
    
    % remove bad samples
    %roiTempChange = nanBadSamples(roiTempChange.',P,Q,zero,1,1).';
    roiTempClean = nanBadSamples(roiTempClean.',P,Q,zero,1,1).';

    allRoiTempClean{s}=roiTempClean;
end

%%

figure;
for s=1:nSubj
    hs(s)=subplot(nSubj,1,s); hold on
    plot(time,mean(allRoiTempClean{s},2),'LineWidth',1.5);
    yl=ylim;
    area([onsetTR*TR offsetTR*TR],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
    ylabel('Temp Change (^oC)');
    xlabel('Time (s)');
    axis tight
    title(titles{s});
end

% % overlay on single axis
% figure; hold on
% for s=1:nSubj
%     hp(s)=plot(time,mean(allRoiTempClean{s},2),'LineWidth',1.5);
%     ylabel('Temp Change (^oC)');
%     xlabel('Time (s)');
%     axis tight
% end
% yl=ylim;
% area([onsetTR*TR offsetTR*TR],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
% htit=title(titlestr);
% set(htit,'FontWeight','normal');
% hlg=legend('sham','active');
% set(hlg,'box','off');


print('-dpng',['../figures/thermo_all_subs' date '.png']);
crop(['../figures/thermo_all_subs' date '.png'],0);



%%
% convert to temperature change


