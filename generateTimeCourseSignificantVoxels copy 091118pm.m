clear all; close all; clc

nTRs=645;
LASERONSETTR=215;
LASEROFFSETTR=430;
testInds=setdiff(1:nTRs,[1:LASERONSETTR  LASEROFFSETTR:nTRs]);
nPerms=1000;
nTRsToSkip=5; % while gradients stabilize
TR=2.8;
time=(0:nTRs-1)*TR;
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
muBoldFilenames={'s8muwBoldEcho1-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    's8muwBoldEcho2-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    's8muwBoldEcho3-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK'};
filenames={'thresh_isSigBefDur_FWHM8_ECHO1-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    'thresh_isSigBefDur_FWHM8_ECHO2-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    'thresh_isSigBefDur_FWHM8_ECHO3-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK'};
brainMaskFilename='resampled_brain_mask+tlrc.BRIK';
lineColor=[0 0 0];
winsize=20; % time series smoothing parameter
legDelRight=0.075;
legDelUp=0.025;
nFilenames=numel(filenames);
laserOrigin = [40 1 34]+[1 1 1]; % LR, AP, IS
% laserOrigin(1)=64-laserOrigin(1);

[~,brainMask,~,~]=BrikLoad(fullfile(pathToData,brainMaskFilename));
brainMask=logical(brainMask);
[NX,NY,NZ]=size(brainMask);
[X,Y,Z]=ndgrid(1:NX,1:NY,1:NZ);
LOx=repmat(laserOrigin(1),[size(X,1) size(X,2) size(X,3)]);
LOy=repmat(laserOrigin(2),[size(X,1) size(X,2) size(X,3)]);
LOz=repmat(laserOrigin(3),[size(X,1) size(X,2) size(X,3)]);
dists=sqrt( (LOx-X).^2 + (LOy-Y).^2 + (LOz-Z).^2 );
% dists1D=dists(:);
% [dists1Dsort,minind]=sort(dists1D,'ascend');
% vox2use=minind(1:214);



hf=figure;
for f=1:nFilenames
    [~,muBold,~,~]=BrikLoad(fullfile(pathToData,muBoldFilenames{f}));
    muBold(:,:,:,1:nTRsToSkip)=0;
    
    [~,isSig,~,~]=BrikLoad(fullfile(pathToData,filenames{f}));
    sigTs=vol2ts(muBold,isSig>0);
    muSigTs=mean(sigTs,2);
    smuSigTs = smoothdata(muSigTs,1,'rloess',winsize);

    % permutation testing
    mockSmuSigTs=zeros(size(muSigTs,1),nPerms);
    for p=1:nPerms
        p
        tmp=surrogateResponseGenerator(muSigTs);
        mockSmuSigTs(:,p)= smoothdata(tmp,1,'rloess',winsize);
    end
    
    [mockSems,mockMus,mockStds]=nansem(mockSmuSigTs,2);
    
    %     nsigTs=vol2ts(muBold,isSig==0 & brainMask & dists<9);
    %     Nsig=size(sigTs,2); %number of significant voxels
    %     Nnsig=size(nsigTs,2); %number of non-significant voxels
    %     vox2use=randperm(Nnsig);
    %     vox2use=vox2use(1:Nsig);
    %     nsigTs=nsigTs(:,vox2use);
    
%     muNsigTs=mean(nsigTs,2);
%     smuNsigTs = smoothdata(muNsigTs,1,'rloess',winsize);
            
    
    hs(f)=subplot(2,3,f);
    hold on

    %hsm2=plot(time(fdr1),smuSigTs1(fdr1),'o','LineStyle','none','MarkerEdgeColor','r','MarkerFaceColor','r');
    %hsm2=plot(time,mockSmuSigTs,'LineWidth',2,'Color','r');
    hsm=plot(time,smuSigTs,'LineWidth',2,'Color',lineColor);
    hmock(f)=shadedErrorBar(time,mockMus,2*mockStds,'k',1);
    set(hmock(f).mainLine,'LineStyle','none');
    set(hmock(f).patch,'FaceColor',[0.7,0.7,0.7]);
    ylim([-0.1 0.1]);
    yl=ylim;
    set(gca,'YTick',[-0.1 0 0.1]);
    harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.25);
    xlabel('Time (mm:ss)','FontSize',20);
    ylabel('BOLD (a.u.)','FontSize',20);
    set(gca,'XTick',[0:300:1800]);
    set(gca,'XTickLabel',{'0:00','5:00','10:00','15:00','20:00','25:00','30:00'});
    xlim([0 1800]);
    axis square
    hlg=legend(harea,'laser on');
    set(hlg,'box','off');
    set(hlg,'FontSize',14);
    %htit(1)=title('Echo 1','FontWeight','normal','FontSize',14);
    
    
end
%
%
% % smuNsigTs1 = smoothdata(muNsigTs1,1,'rloess',winsize);
%
% [~,muBold2,~,~]=BrikLoad(muBoldFilename2);
% [~,muBold3,~,~]=BrikLoad(muBoldFilename3);
%
% % throw out first 5 TRs
% muBold1(:,:,:,1:nTRsToSkip)=0;
% muBold2(:,:,:,1:nTRsToSkip)=0;
% muBold3(:,:,:,1:nTRsToSkip)=0;
%
% [~,isSig1,~,~]=BrikLoad(filenames{1});
% [~,isSig2,~,~]=BrikLoad(filenames{2});
% [~,isSig3,~,~]=BrikLoad(filenames{3});
%
% sigTs1=vol2ts(muBold1,isSig1>0);
% sigTs2=vol2ts(muBold2,isSig2>0);
% sigTs3=vol2ts(muBold3,isSig3>0);
%
% % nsigTs1=vol2ts(muBold1,isSig1==0);
% % nsigTs2=vol2ts(muBold2,isSig2==0);
%
% muSigTs1=mean(sigTs1,2);
% muSigTs2=mean(sigTs2,2);
% muSigTs3=mean(sigTs3,2);
% %
% % muNsigTs1=mean(nsigTs1,2);
% % muNsigTs2=mean(nsigTs2,2);
%
% %% smooth data
% smuSigTs1 = smoothdata(muSigTs1,1,'rloess',winsize);
% smuSigTs2 = smoothdata(muSigTs2,1,'rloess',winsize);
% smuSigTs3 = smoothdata(muSigTs3,1,'rloess',winsize);
%
% % smuNsigTs1 = smoothdata(muNsigTs1,1,'rloess',winsize);
% % smuNsigTs2 = smoothdata(muNsigTs2,1,'rloess',winsize);
% % smuNsigTs3 = smoothdata(muNsigTs3,1,'rloess',winsize);
%
% %% permutation testing
% nsmuTs1=zeros(size(smuSigTs1,1),nPerms);
% nsmuTs2=zeros(size(smuSigTs2,1),nPerms);
% nsmuTs3=zeros(size(smuSigTs3,1),nPerms);
% for p=1:nPerms
%     p
%     tmp1=surrogateResponseGenerator(muSigTs1);
%     tmp2=surrogateResponseGenerator(muSigTs2);
%     tmp3=surrogateResponseGenerator(muSigTs3);
%     nsmuTs1(:,p)= smoothdata(tmp1,1,'rloess',winsize);
%     nsmuTs2(:,p) = smoothdata(tmp2,1,'rloess',winsize);
%     nsmuTs3(:,p) = smoothdata(tmp3,1,'rloess',winsize);
% end
%
% %%  compute p-values
% pvals1=mean(nsmuTs1>repmat(smuSigTs1,[1 nPerms]),2);
% pvals2=mean(nsmuTs2>repmat(smuSigTs2,[1 nPerms]),2);
% pvals3=mean(nsmuTs3>repmat(smuSigTs3,[1 nPerms]),2);
%
% [t1,tfdr1]=fdr(pvals1(testInds),0.05);
% [t2,tfdr2]=fdr(pvals2(testInds),0.05);
% [t3,tfdr3]=fdr(pvals3(testInds),0.05);
%
% fdr1=false(nTRs,1);  fdr1(testInds)=tfdr1==1;
% fdr2=false(nTRs,1);  fdr2(testInds)=tfdr2==1;
% fdr3=false(nTRs,1);  fdr3(testInds)=tfdr3==1;
% %%





sublabel(hs,50,-20,'FontSize',16,'FontWeight','Bold');
print -dpng ../figures/timeSeriesFigure.png
crop('../figures/timeSeriesFigure.png',0);