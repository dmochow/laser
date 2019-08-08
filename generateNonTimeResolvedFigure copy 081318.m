% 06/04/18
% Answer the basic question: is there an effect
clear all; close all; clc

%% the precomputed data
dataFilenames={'../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-TSHIFT1-Echo1-05-Jun-2018'; ...
    '../data/precomputed/allControlBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-TSHIFT1-Echo1-04-Jun-2018'};
nFilenames=numel(dataFilenames);
figFilename='nonTimeResolvedFigure.png';
% dataFilenames={'../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo1-02-Jun-2018.mat'; ...
%     '../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo2-31-May-2018.mat'; ...
%     '../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo3-31-May-2018.mat'};

%% constants
TR=2.8; nTR=645; time=(0:nTR-1)*TR;
laserOnsetTR=215;
laserOffsetTR=430;
nlags=5; % optimized in optimizeARmodelOrder.m
tPre=[0 600];
tLaser=[600 1200];
tPost=[1200 1800];
alpha=0.05; % sig level
del=0.01; % sig bar height
delTy=0.05;  % sig text up
delTx=50; % sig text left
delLgUp=0.05; % legend up
delLgRt=0.12; % legend right
delUp=0.05; % subplots 2 and 3 up
delBoldUp=0.075; % bring time series plots up to reduce white space
%% load data
figure;
for e=1:nFilenames
    load(dataFilenames{e},'alloBolds');
    allBolds2D=cat(2,alloBolds{:});
    [semBold,grandBold]=nansem(allBolds2D,2);
    
    y=grandBold;
    X=tplitz(y,nlags); X=X(:,2:end-1);
    
    % is pre different from peri?
    [h1,p1] = chowtest(X(1:laserOffsetTR,:),y(1:laserOffsetTR),laserOnsetTR);
    
    % is pre different from post?
    [h2,p2] = chowtest(cat(1,X(1:laserOnsetTR,:),X(laserOffsetTR+1:end,:)),cat(1,y(1:laserOnsetTR),y(laserOffsetTR+1:end)),laserOnsetTR);
    
    % draw
    hs(e)=subplot(nFilenames,2,2+e); hold on
    hsh=shadedErrorBar(time,grandBold,semBold); yl=ylim;
    
    harea0=area(tPre,yl(2)*[1 1],'BaseValue',yl(1));
    set(harea0,'FaceColor',[0.7 0.7 0.7]);
    set(harea0,'FaceAlpha',0.25);
    set(harea0,'EdgeColor','none');
    
    harea1=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
    set(harea1,'FaceColor',[0.35 0.35 0.35]);
    set(harea1,'FaceAlpha',0.25);
    set(harea1,'EdgeColor','none');
    
    harea2=area(tPost,yl(2)*[1 1],'BaseValue',yl(1));
    set(harea2,'FaceColor',[0.7 0.7 0.7]);
    set(harea2,'FaceAlpha',0.25);
    set(harea2,'EdgeColor','none');
    
    
    lineHt=yl(2)+del;
    lineHt2=(yl(2))*1.5+del;
    if p1<alpha
        plot([300 300],[lineHt lineHt+del],'k');
        plot([300 900],[lineHt+del lineHt+del],'k');
        plot([900 900],[lineHt lineHt+del],'k');
        text(600-delTx,lineHt+delTy,sprintf(['p=%0.3f'],p1));
    else
        plot([300 300],[lineHt lineHt+del],'k');
        plot([300 900],[lineHt+del lineHt+del],'k');
        plot([900 900],[lineHt lineHt+del],'k');
        text(600-delTx,lineHt+delTy,'n.s.');
    end
    
    if p2<alpha
        plot([300 300],[lineHt2 lineHt2+del],'k');
        plot([300 1500],[lineHt2+del lineHt2+del],'k');
        plot([1500 1500],[lineHt2 lineHt2+del],'k');
        text(900-delTx,lineHt2+delTy,sprintf(['p=%0.3f'],p2));
    else
        plot([300 300],[lineHt2 lineHt2+del],'k');
        plot([300 1500],[lineHt2+del lineHt2+del],'k');
        plot([1500 1500],[lineHt2 lineHt2+del],'k');
        text(900-delTx,lineHt2+delTy,'n.s.');
    end
    
    if e==1, hlg=legend([harea1 harea2],'Laser On','Laser Off'); end
end

set(hs(:),'XLim',[0 1801]);
set(hs(:),'Xtick',[0 600 1200 1800]);
set(hs(1),'XTickLabel',{'0:00','10:00','20:00','30:00'});
set(hs(2),'XTickLabel',{'0:00','10:00','20:00','30:00'});
set(get(hs(1),'xlabel'),'String','Time (mm:ss)');
set(get(hs(2),'xlabel'),'String','Time (mm:ss)');

%set(hs(:),'Ytick',[-0.2 0 0.2]);
set(hs(:),'Ytick',[-0.2 -0.1 0 0.1 0.2]);
set(hs(2),'YtickLabel',{'','','','',''});
set(get(hs(1),'ylabel'),'String','BOLD (a.u.)');
set(get(hs(2),'ylabel'),'String','');

set(get(hs(1),'Title'),'String','ROI');
set(get(hs(2),'Title'),'String','Contralateral');
set(get(hs(1),'Title'),'FontWeight','Bold');
set(get(hs(2),'Title'),'FontWeight','Bold');


lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)+delLgRt lgPos(2)+delLgUp lgPos(3) lgPos(4)]);
set(hlg,'Box','off');

%%
art=imread('../figures/mrilaser.jpg');
hss(4)=subplot(221);
imagesc(art); axis off


%% draw anatomy with overlaid ROI
[~,anat,~,~]=BrikLoad('../data/S16/NII/anat+orig');
%nii=load_untouch_nii('../data/S16/NII/anat.nii'); anat2=nii.img;
[~,roi,~,~]=BrikLoad('../data/S16/NII/roi+orig');
axsl=146; sagsl=137; corsl=256-56;

% im1=squeeze(anat(:,:,axsl)); 
% im2=squeeze(roi(:,:,axsl));
% im1=cat(3,im1,im1,im1);
% im2=cat(3,im2,im2,im2);
% axImg = imfuse(im1,im2,'blend');

axImg = rgb2gray(imfuse(squeeze(anat(:,:,axsl)),squeeze(roi(:,:,axsl))));
sagImg = rgb2gray(imfuse(squeeze(anat(sagsl,:,:)),squeeze(roi(sagsl,:,:))));
corImg = rgb2gray(imfuse(squeeze(anat(:,corsl,:)),squeeze(roi(:,corsl,:))));

hss(1)=subplot(264);
imagesc(axImg.'); axis square; axis xy;  colormap bone; axis off;
hss(2)=subplot(265);
imagesc(sagImg.'); axis square; axis xy; colormap bone; axis off; 
hss(3)=subplot(266);
imagesc(corImg.'); axis square; axis xy;  colormap bone; axis off;

% resize mri images
scFc=1.5;
mDelDwn=0.08;
mDelLeft=0.025;
mDelRight=0.025;
sagRange=[30, 224-30];
corRange=[30, 256-30];
axRange=[30, 256-30];
%%
for i=1:3
    mpos=get(hss(i),'Position');
    
    switch i
        case 1
            set(hss(i),'Position',[mpos(1)-mDelLeft mpos(2)-mDelDwn scFc*mpos(3) scFc*mpos(4)]);
            set(hss(i),'XLim',sagRange);
            set(hss(i),'YLim',corRange);
        case 2
             set(hss(i),'Position',[mpos(1) mpos(2)-mDelDwn scFc*mpos(3) scFc*mpos(4)]);
            set(hss(i),'XLim',corRange);
            set(hss(i),'YLim',axRange);
        case 3
            set(hss(i),'Position',[mpos(1)+mDelRight mpos(2)-mDelDwn scFc*mpos(3) scFc*mpos(4)]);
            set(hss(i),'XLim',sagRange);
            set(hss(i),'YLim',axRange);
    end
end

%% time series plot up
for i=1:2
    spos=get(hs(i),'Position');
    set(hs(i),'Position',[spos(1) spos(2)+delBoldUp spos(3) spos(4)]);
end

%%
htxt(1)=text(-775,275,'A','FontSize',20);
htxt(2)=text(-325,275,'B','FontSize',20);
htxt(3)=text(-775,-80,'C','FontSize',20);
htxt(4)=text(-325,-80,'D','FontSize',20);


%% print to file
print('-dpng',figFilename);% ../figures/nonTimeResolvedFigure
crop(figFilename);





% %% Analysis 3: split at 2-minute windows following laser onset
% %% set up model and data structures
% nlags=5;
% bps2=TRonset;
% bpsTime2=bps2*TR;
% wins=[(10:1:29)' (11:1:30)'];
% wins=round(wins*60/TR);
% winsTime=wins*TR;
% nWins=size(wins,1);
% pvals2=zeros(nWins,1);
% for w=1:nWins
%     w
%     y=cat(1,grandBold(1:TRonset),grandBold(wins(w,1):wins(w,2)));
%     y=surrogateResponseGenerator(y);
%     X=tplitz(y,nlags); X=X(:,2:end-1);
%     [h,p] = chowtest(X,y,bps2);
%     pvals2(w)=p;
% end
% % run Chow tests
%
% isSig2=pvals2<alpha;
% % [pfdr2,isSig2]=fdr(pvals2,alpha);
% % notSig=~isSig;
%
% % draw figure
% figure;
% hs(1)=subplot(211);
% hsh=shadedErrorBar(time,grandBold,semBold);
% hs(2)=subplot(212);
% hp=plot(winsTime(:,1),pvals2,'*k'); hold on
% plot(winsTime(isSig2),pvals2(isSig2),'*b','MarkerFaceColor','b');
% %plot(bpsTime(notSig),pvals(notSig),'or','MarkerFaceColor','r');
% ylabel('P-value');
% xlabel('Time (s)');
% set(hs(1),'Xtick',[0 600 1200 1800]);
% set(hs(2),'Xtick',[0 600 1200 1800]);
% set(hs(2),'Xlim',[0 1800]); set(hs(2),'Ylim',[0 1]);
% % Mdl = fitlm(X,y);
% %
% % %
% % res = Mdl.Residuals.Raw;
% % figure;
% % subplot(221);
% % plotResiduals(Mdl,'lagged');
% % subplot(222);
% % plotResiduals(Mdl,'caseorder');
%
