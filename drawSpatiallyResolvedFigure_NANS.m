% 08/07/18
% draw some colored brains
clear all; close all; clc

%%
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
brainMaskFilename=fullfile(pathToData,'brain_mask+tlrc.BRIK');
roiMaskFilename=fullfile(pathToData,'mu_roi_r19_z26+tlrc.BRIK');
cropMri=1; % whether to crop each colored mri plot
widthMultiplier=1.31;
heightMultiplier=1.31;
figw=8.5;
figh=8.5;
subx_start=0.025;
subx_end=0.975;
subx_total_width=subx_end-subx_start;
for e=1:3
    isSigBefDurFilename=fullfile(pathToData,['resampled_thresh_isSigBefDur_FWHM8_ECHO' num2str(e) '-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK']);
    isSigBefAftFilename=fullfile(pathToData,['resampled_thresh_isSigBefAft_FWHM8_ECHO' num2str(e) '-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK']);
    
    %% load the anatomy and the significance clusters
    [~, thisIsSigBefDur, ~, ~] = BrikLoad (isSigBefDurFilename);
    [~, thisIsSigBefAft, ~, ~] = BrikLoad (isSigBefAftFilename);
    
    %% reflect for pcolor's sake
    isSigBefDur(:,:,:,e)=thisIsSigBefDur(end:-1:1,:,:);
    isSigBefAft(:,:,:,e)=thisIsSigBefAft(end:-1:1,:,:);
end
[~, anat, ~, ~] = BrikLoad (anatFilename);
anat=anat(end:-1:1,:,:);

[~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
brainMask=brainMask(end:-1:1,:,:);

[~, roiMask, ~, ~] = BrikLoad (roiMaskFilename);
roiMask=roiMask(end:-1:1,:,:);
roiMask=logical(roiMask);

%%
nEchos=3;
nRows=nEchos*2; % for each echo, before vs during, before vs after
axialSliceInds=78:2:96;
nCols=numel(axialSliceInds);
echoInds=repmat((1:nEchos)',2,1);
subx_width=subx_total_width/nCols-0.02;

hf=figure('Units','inches','Position',[0 0 figw figh]);
for c=1:nCols
    for e=1:nEchos
        
        thisIsSigBefDur=isSigBefDur(:,:,:,e);
        thisIsSigBefAft=isSigBefAft(:,:,:,e);
        
        hs(e,c)=subplot(nRows,nCols,(e-1)*nCols+c);
        
        % 
        %h = colorBrain(anat(:,:,axialSliceInds(c)),roiMask(:,:,axialSliceInds(c)),thisIsSigBefDur(:,:,axialSliceInds(c))==1,brainMask(:,:,axialSliceInds(c)),cropMri);
        
        % no activation shown
        h = colorBrain(anat(:,:,axialSliceInds(c)),roiMask(:,:,axialSliceInds(c)),thisIsSigBefDur(:,:,axialSliceInds(c))==2,brainMask(:,:,axialSliceInds(c)),cropMri);
        
        pos1=get(hs(e,c),'Position');
        aspRatio=pos1(4)/pos1(3);
        set(hs(e,c),'Position',[subx_start+subx_width*(c-1) pos1(2) subx_width subx_width*aspRatio]);
        
        if e==1 % no sig voxels at echos 2 and 3
            hs(e+nEchos,c)=subplot(nRows,nCols,(e+nEchos-1)*nCols+c);
            
            %h = colorBrain(anat(:,:,axialSliceInds(c)),roiMask(:,:,axialSliceInds(c)),thisIsSigBefAft(:,:,axialSliceInds(c))==1,brainMask(:,:,axialSliceInds(c)),cropMri);
            
            h = colorBrain(anat(:,:,axialSliceInds(c)),roiMask(:,:,axialSliceInds(c)),thisIsSigBefAft(:,:,axialSliceInds(c))==2,brainMask(:,:,axialSliceInds(c)),cropMri);
           
            
            pos1=get(hs(e+nEchos,c),'Position');
            set(hs(e+nEchos,c),'Position',[subx_start+subx_width*(c-1) pos1(2) subx_width subx_width*aspRatio]);
            htit(e,c)=title(sprintf(['z=' num2str(axialSliceInds(c))]),'FontWeight','normal');
        end
        
    end
end

% move all rows down and to the right
delDown=0.1; delRight=0.05;
for row=1:1+nEchos
    for col=1:nCols
        pos=get(hs(row,col),'Position');
        set(hs(row,col),'Position',[pos(1)+delRight pos(2)-delDown pos(3) pos(4)]);
    end
end

% move fourth row down even more
for col=1:nCols
    pos=get(hs(1+nEchos,col),'Position');
    set(hs(1+nEchos,col),'Position',[pos(1) pos(2)-delDown pos(3) pos(4)]);
end

%%
% add timing figures
htim(1)=subplot(9,9,73);yl=ylim; xlim=[0 3];
ylim([0 3]);
hold on
harea=area([2 3],yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0 0 0]);
set(harea,'FaceAlpha',0.25);
set(harea,'EdgeColor','none');

harea=area([1 2],yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.25);
set(harea,'EdgeColor','none');

harea=area([0 1],yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0 0 0]);
set(harea,'FaceAlpha',0.25);
set(harea,'EdgeColor','none');
text(0.25,0.85,'off');
text(1.25,0.85,'on');
text(2.25,0.85,'off');

%set(gca,'XTick',[0 3]);
%set(gca,'XTickLabel',{'0:00','30:00'});

axis off; axis tight;

plot([0.5 0.5],[1.1 1.2],'k')
plot([0.5 1.5],[1.2 1.2],'k')
plot([1.5 1.5],[1.1 1.2],'k')
text(0.8,1.45,'vs')
pos1=get(htim(1,1),'Position');
pos2=get(hs(1,1),'Position');
set(htim(1),'Position',[pos2(1)+0.05 pos2(2)+0.115 pos1(3)*2 pos1(4)*1]);

xx=6.8; delyy=2.6;
text(xx,0.25,'Echo 1');
text(xx,0.25-delyy,'Echo 2');
text(xx,0.25-2*delyy,'Echo 3');
text(xx,-9.2,'Echo 1');
%%
% second timing figure
htim(2)=subplot(9,9,74);yl=ylim; xlim=[0 3];
ylim([0 3]);
hold on
harea=area([2 3],yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0 0 0]);
set(harea,'FaceAlpha',0.25);
set(harea,'EdgeColor','none');

harea=area([1 2],yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.25);
set(harea,'EdgeColor','none');

harea=area([0 1],yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0 0 0]);
set(harea,'FaceAlpha',0.25);
set(harea,'EdgeColor','none');
text(0.25,0.85,'off');
text(1.25,0.85,'on');
text(2.25,0.85,'off');

axis off; axis tight;

plot([0.5 0.5],[1.1 1.2],'k')
plot([0.5 2.5],[1.2 1.2],'k')
plot([2.5 2.5],[1.1 1.2],'k')
text(1.35,1.4,'vs')
pos1=get(htim(2),'Position');
pos2=get(hs(4,1),'Position');
set(htim(2),'Position',[pos2(1)+0.05 pos2(2)+0.115 pos1(3)*2 pos1(4)*1]);

%%
% add z labels to final row
for c=1:nCols, moveTitle(htit(1,c),0,250,0); end

%%
% add mock patches for legend
x = [0 eps eps 0];
y = [0 0 eps eps];
hpatch1=patch(x,y,'g');
hpatch2=patch(x,y,'r');
set(hpatch1,'EdgeColor','none');
set(hpatch2,'EdgeColor','none');
hlg=legend([hpatch1],'Light absorption');
%hlg=legend([hpatch1 hpatch2],'Light absorption','Significant BOLD change');
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)+0.525 lgPos(2)+0.565 lgPos(3) lgPos(4)]);
set(hlg,'FontSize',16);
set(hlg,'Box','off');


%%
% sublabel
sublabel([hs(1,1) hs(4,1)],-30,0,'FontSize',16,'FontWeight','Bold');
print -dpng ../figures/figure3a
crop('../figures/figure3a.png');

%%
% Jaccard
nEchos=3;
s1=roiMask;     
inds1=find(s1);
J1=zeros(nEchos,1); 
J2=zeros(nEchos,1);
for e=1:nEchos
    s2=(squeeze(isSigBefDur(:,:,:,e)))>0;
    inds2=find(s2);
    %J1(e)=sum(intersect(inds1,inds2))/sum(union(inds1,inds2));
    J1(e)=sum(intersect(inds1,inds2))/sum(inds2); % % sig voxels that are in the ROI
    
    s2=(squeeze(isSigBefAft(:,:,:,e)))>0;
    inds2=find(s2);
    %J2(e)=sum(intersect(inds1,inds2))/sum(union(inds1,inds2));
    J2(e)=sum(intersect(inds1,inds2))/sum(inds2);
end
J1
J2
