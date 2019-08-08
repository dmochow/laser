clear all; close all; clc
% read in filtered data from FSL and try to analyze it
% after "denoising" with fslregfilt
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));
addpath(genpath('./vol3d'));

%%
subjIndx=1;
alpha=0.05; % for output nii
nPerms=500; % for significance testing
%%
switch subjIndx
    case 1
        maskFilename='../data/S01/FSL.ica/filtered_func_data.ica/mask.nii';
        anatFilename='../data/S01/FSL.orig/T1w_MPR_BIC_v1.nii';
        dataFilename='../data/S01/FSL.orig/BOLD_AP.ica/denoised_data.nii';
        tOn=401:800;
        tOff=[1:400 801:1200];
        tOff_1=1:400;
        tOff_2=801:1200;
        fs=1/1.5;  % 
%         outFilename='../data/S01/FSL.ica/signalRatio_whighpass.nii';
        outFilename2='../data/S01/FSL.orig/BOLD_AP.ica/design_pvals_updown.nii';;
    case 2
        maskFilename='../data/S02/FSL.ica/filtered_func_data.ica/mask.nii';
        dataFilename='../data/S02/FSL.ica/filtered_func_data.nii';
        tOn=216:430;
        tOff=[1:215 431:645];
        fs=1/2.8;  
        outFilename='../data/S02/FSL.ica/signalRatio_whighpass.nii';
        outFilename2='../data/S02/FSL.ica/pvalues_whighpass.nii';
end

%%
maskNii=load_untouch_nii(maskFilename);
mask=logical(maskNii.img);

%%
dataNii=load_untouch_nii(dataFilename);
data=double(dataNii.img);
data_4D=data;
data=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
brainData=data(mask(:),:);

%%
% for each voxel in the brain, compute fit of time series to the design
% function (flat-ramp-flat)
%%
design=zeros(numel(tOff_1)+numel(tOn)+numel(tOff_2),1);
design(tOn)=linspace(0,1,numel(tOn));
design(tOff_2)=linspace(1,0,numel(tOff_2));
nBrainVoxels=numel(find(mask(:)));
nPerm=500;
r2=zeros(nBrainVoxels,1);
r2_mock=zeros(nBrainVoxels,nPerm);
tic
for v=1:nBrainVoxels
    v
    Y=brainData(v,:)';
    X=design;
    [b,bint,r,rint,stats] = regress(Y,[X ones(size(X,1),1)]);
    r2(v)=stats(1);
    
    for p=1:nPerm
        Y_mock=surrogateResponseGenerator(Y);
        [b,bint,r,rint,stats_mock] = regress(Y_mock,[X ones(size(X,1),1)]);
        r2_mock(v,p)=stats_mock(1);
        
    end
    
end
toc



% 
% %%
% % basic figure
% clrs=[0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250];

% 
% 
% % %%
% % % permutation testing
% % mock_signalRatio=zeros(size(signalRatio,1),nPerms);
% % nOn=numel(tOn);
% % nOff=numel(tOff);
% % for p=1:nPerms
% %     p
% %     
% %     mock_brainData=surrogateResponseGenerator(brainData')';
% %     
% %     mock_signalOn=std(mock_brainData(:,tOn),[],2);
% %     mock_signalOff=std(brainData(:,tOff),[],2);
% %     mock_signalRatio(:,p)=mock_signalOn./mock_signalOff;
% %     
% % end
% % 
% % %% 
% compute pvals
%nVoxels=size(signalRatio,1);
%%
pvals=zeros(nBrainVoxels,1);
for v=1:nBrainVoxels
    pvals(v)=mean(r2_mock(v,:)>r2(v));
end

% correct p-values for multiple comparisons
[p_fdr, p_masked] = fdr( pvals, alpha);
%%

% % 
% % 
% % %%
%%
% tmp=zeros(size(mask(:)));
% tmp(mask(:))=r2;
% r2Volume=reshape(tmp,[size(mask,1),size(mask,2),size(mask,3)]);
% figure;
% subplot(121); vol3d('cdata',r2Volume);
%%


% % 
tmp=zeros(size(mask(:)));
tmp(mask(:))=p_masked;
pvalVolume=reshape(tmp,[size(mask,1),size(mask,2),size(mask,3)]);
% % 
% % 
% % 
% % % ratioVolume_withNan=ratioVolume;
% % % ratioVolume_withNan(~mask)=NaN;
% % % ratiosLeft=ratioVolume_withNan(1:45,:,:); nanmean(ratiosLeft(:))
% % % ratiosRight=ratioVolume_withNan(46:end,:,:);nanmean(ratiosRight(:))
% % % ratiosFront=ratioVolume_withNan(:,1:45,:); nanmean(ratiosFront(:))
% % % ratiosBack=ratioVolume_withNan(:,46:end,:);nanmean(ratiosBack(:))
% % % ratiosUp=ratioVolume_withNan(:,:,31:60); nanmean(ratiosUp(:))
% % % ratiosDown=ratioVolume_withNan(:,:,1:30);nanmean(ratiosDown(:))
% % 
% % %%

%subplot(122); vol3d('cdata',pvalVolume);
% axis tight;  daspect([1 1 .4])
% alphamap('rampup');
% alphamap(.5 .* alphamap);
% 
% % %%
% % % write output nii
% % outNii=maskNii;
% % outNii.img=ratioVolume;
% % save_untouch_nii(outNii,outFilename);
% % 
%%
outNii2=maskNii;
outNii2.img=pvalVolume;
save_untouch_nii(outNii2,outFilename2);

%%
% 
% 
% %%
% % figure;
% % plot(mean(brainData.^2,1))
% 
% %%
% vind=[36,73,22]+1;
% vx_right=33:43; % right
% vx_left=47:57; % left
% vy_front=63:73;
% vy_back=11:21;
% vz=28:38;
% time=(1:1200)*1.5;
% 
% subdata_left=data_4D(vx_left,vy_front,vz,:);
% subdata_right=data_4D(vx_right,vy_front,vz,:);
% subdata_back=data_4D(vx_right,vy_back,vz,:);
% subdata_cross=data_4D(vx_left,vy_back,vz,:);
% ts_left=reshape(subdata_left,[numel(vx)*numel(vy)*numel(vz), size(subdata,4)]);
% ts_right=reshape(subdata_right,[numel(vx)*numel(vy)*numel(vz), size(subdata,4)]);
% ts_back=reshape(subdata_back,[numel(vx)*numel(vy)*numel(vz), size(subdata,4)]);
% ts_cross=reshape(subdata_cross,[numel(vx)*numel(vy)*numel(vz), size(subdata,4)]);
% 
% width=600;
% 
% 
% 
% figure;
% subplot(221);
% plot(time,mean(ts_right)); title('right front (target)')
% xlim([time(1) time(end)]);
% lim=axis;
% x=[lim(1) lim(1) 600 600];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOn=patch(x,y,[.8 .8 .8]);
% set(hOn,'FaceAlpha',0.3,'LineStyle','none');
% x=[1200 1200 1800 1800];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOff=patch(x,y,[.8 .8 .8]);
% set(hOff,'FaceAlpha',0.3,'LineStyle','none');
% htxtOff=text(250,[lim(4)-lim(3)]*0.75+lim(3),'off');
% htxtOn=text(850,[lim(4)-lim(3)]*0.25+lim(3),'on');
% htxtOff2=text(1450,[lim(4)-lim(3)]*0.75+lim(3),'off');
% xlabel('Time (s)');
% ylabel('BOLD (a.u.)');
% 
% subplot(222);
% plot(time,mean(ts_left)); title('left front');
% xlim([time(1) time(end)]);
% lim=axis;
% x=[lim(1) lim(1) 600 600];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOn=patch(x,y,[.8 .8 .8]);
% set(hOn,'FaceAlpha',0.3,'LineStyle','none');
% x=[1200 1200 1800 1800];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOff=patch(x,y,[.8 .8 .8]);
% set(hOff,'FaceAlpha',0.3,'LineStyle','none');
% 
% 
% 
% 
% subplot(223);
% plot(time,mean(ts_back)); title('right back');
% xlim([time(1) time(end)]);
% lim=axis;
% x=[lim(1) lim(1) 600 600];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOn=patch(x,y,[.8 .8 .8]);
% set(hOn,'FaceAlpha',0.3,'LineStyle','none');
% x=[1200 1200 1800 1800];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOff=patch(x,y,[.8 .8 .8]);
% set(hOff,'FaceAlpha',0.3,'LineStyle','none');
% 
% 
% 
% subplot(224);
% plot(time,mean(ts_cross)); title('left back');
% xlim([time(1) time(end)]);
% lim=axis;
% x=[lim(1) lim(1) 600 600];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOn=patch(x,y,[.8 .8 .8]);
% set(hOn,'FaceAlpha',0.3,'LineStyle','none');
% x=[1200 1200 1800 1800];
% y=[lim(3) lim(4) lim(4) lim(3)];
% hOff=patch(x,y,[.8 .8 .8]);
% set(hOff,'FaceAlpha',0.3,'LineStyle','none');
% 
% print -dpng ../figures/mean_signal_4ROI
% 
% %%
% 
% anatNii=load_untouch_nii(anatFilename);
% anat=double(anatNii.img);
% figure;
% subplot(221)
% imagesc([0 224],[0 256],anat(:,end:-1:1,155).'); hold on
% x=[70 70 95 95];
% y=[80 55 55 80];
% h=patch(x,y,[0.9 0 0],'FaceAlpha',0.5,'LineStyle','none');
% colormap(bone);
% title('ROI: right front');
% axis off
% 
% subplot(222)
% imagesc([0 224],[0 256],anat(:,end:-1:1,155).'); hold on
% x=[115 115 140 140];
% y=[80 55 55 80];
% h=patch(x,y,[0.9 0 0],'FaceAlpha',0.5,'LineStyle','none');
% colormap(bone);
% title('ROI: left front');
% axis off
% 
% subplot(223)
% imagesc([0 224],[0 256],anat(:,end:-1:1,155).'); hold on
% x=[70 70 95 95];
% y=[195 220 220 195];
% h=patch(x,y,[0.9 0 0],'FaceAlpha',0.5,'LineStyle','none');
% colormap(bone);
% title('ROI: right back');
% axis off
% 
% subplot(224)
% imagesc([0 224],[0 256],anat(:,end:-1:1,155).'); hold on
% x=[115 115 140 140];
% y=[195 220 220 195];
% h=patch(x,y,[0.9 0 0],'FaceAlpha',0.5,'LineStyle','none');
% colormap(bone);
% title('ROI: left back');
% axis off
% 
% print -dpng ../figures/ROIs
% 
% %%
% % subdata_right_on=subdata_right(:,:,:,tOn);
% % mu_right_on=mean(subdata_right_on(:));
% % sem_right_on=std(subdata_right_on(:)); %/sqrt(numel(subdata_right_on(:)));
% % 
% % subdata_right_off_1=subdata_right(:,:,:,tOff_1);
% % mu_right_off_1=mean(subdata_right_off_1(:));
% % sem_right_off_1=std(subdata_right_off_1(:)); %/sqrt(numel(subdata_right_on(:)));
% % 
% figure
% 
% subplot(221)
% 
% Y=mean(ts_right(:,tOn),1).';
% X=(tOn(1):1:tOn(end)).';
% B_on = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_right(:,tOff_1),1).';
% X=(tOff_1(1):1:tOff_1(end)).';
% B_off_1 = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_right(:,tOff_2),1).';
% X=(tOff_2(1):1:tOff_2(end)).';
% B_off_2 = regress(Y,[X ones(size(Y,1),1)]);
% 
% hbar(1)=bar([1 2 3],[B_off_1(1) B_on(1) B_off_2(1)]);
% ylim([-0.15 0.15])
% title('right front (target)');
% ylabel('slope');
% set(gca,'XTickLabel',{'off','on','off'})
% 
% subplot(222)
% 
% Y=mean(ts_left(:,tOn),1).';
% X=(tOn(1):1:tOn(end)).';
% B_on = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_left(:,tOff_1),1).';
% X=(tOff_1(1):1:tOff_1(end)).';
% B_off_1 = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_left(:,tOff_2),1).';
% X=(tOff_2(1):1:tOff_2(end)).';
% B_off_2 = regress(Y,[X ones(size(Y,1),1)]);
% 
% hbar(2)=bar([1 2 3],[B_off_1(1) B_on(1) B_off_2(1)]);
% ylim([-0.15 0.15])
% title('left front');
% set(gca,'XTickLabel',{'off','on','off'})
% 
% subplot(223)
% 
% Y=mean(ts_back(:,tOn),1).';
% X=(tOn(1):1:tOn(end)).';
% B_on = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_back(:,tOff_1),1).';
% X=(tOff_1(1):1:tOff_1(end)).';
% B_off_1 = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_back(:,tOff_2),1).';
% X=(tOff_2(1):1:tOff_2(end)).';
% B_off_2 = regress(Y,[X ones(size(Y,1),1)]);
% 
% hbar(3)=bar([1 2 3],[B_off_1(1) B_on(1) B_off_2(1)]);
% ylim([-0.15 0.15])
% title('right back');
% set(gca,'XTickLabel',{'off','on','off'})
% 
% subplot(224)
% 
% Y=mean(ts_cross(:,tOn),1).';
% X=(tOn(1):1:tOn(end)).';
% B_on = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_cross(:,tOff_1),1).';
% X=(tOff_1(1):1:tOff_1(end)).';
% B_off_1 = regress(Y,[X ones(size(Y,1),1)]);
% 
% Y=mean(ts_cross(:,tOff_2),1).';
% X=(tOff_2(1):1:tOff_2(end)).';
% B_off_2 = regress(Y,[X ones(size(Y,1),1)]);
% 
% hbar(4)=bar([1 2 3],[B_off_1(1) B_on(1) B_off_2(1)]);
% ylim([-0.15 0.15])
% title('left back');
% set(gca,'XTickLabel',{'off','on','off'})
% 
% 
% set(hbar(:),'FaceColor',[0.8 0.8 0.8]);
% set(hbar(:),'LineStyle','none');
% 
% print -dpng ../figures/slopes_4ROI
