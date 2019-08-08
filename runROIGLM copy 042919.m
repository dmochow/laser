clear all; close all; clc

%% TODO: bring in laser timing from biopac
%% TODO: add second regression for post-laser

PRINT_FIGURE=0;

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

subjectNames={'Greg','Emircan','Pat','Flores','Christian',...
    'Jamal','Kevin','Destiny','Jacek','Marcelo',...
'James','Enisena','Felicia','Ridmila','Karem',...
'Blerita','Duc','Anu','Javier','Heninstoa'};

% the three time segments
t1=[1:215];
t2=[216:430];
t3=[431:645];
X=zeros(size([t1 t2 t3]));
X(t2)=1;
X=X(:);
Xnull=ones(size(X));
alpha=0.05;
gmIndex=2;


for e=1:3
    echoStr=num2str(e);
    filename=['shs8oBoldEcho' echoStr '+tlrc.BRIK'];
    %filename=['pscBoldEcho' echoStr '+tlrc.BRIK']
    
    for s=1:nSubjects
        [e,s]
        pathToData=['../data/' subjStrs{s} '/NII/']
        
        [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,filename));
        
%         brainMaskFilename='resampled_brain_mask+tlrc';
%         [~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
%         brainMask=logical(brainMask);
        
        % grey matter mask
        resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
        [~, seg, info, ~] = BrikLoad (resampledSegFilename);
        brainMask=logical(seg==gmIndex);
        
        roiMaskFilename='resampled_roi_r21_z39+tlrc.BRIK';
        [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
        roiMask=logical(roiMask);
        
        if isempty(roiMask) % not all subjects had markers
            roiMaskFilename='resampled_mu_roi_r21_z39+tlrc.BRIK';
            [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
            roiMask=logical(roiMask);
        end
        
        finalMask=roiMask&brainMask;
        
        %%
        ts=vol2ts(bold,finalMask);
        allts{e,s}=ts;
        nVoxels=size(ts,2);
        for v=1:nVoxels
            statsf=myGLM(ts(:,v),X,1); % full model
            statsr=myGLM(ts(:,v),Xnull,0); % reduced model
            allBetas{e,s}(v)=statsf.B(1);
            ssef=statsf.sse;
            doff=statsf.dof;
            sser=statsr.sse;
            dofr=statsr.dof;
            Fstat=((sser-ssef)/(dofr-doff)) / (ssef/doff);
            pval{s}(v)=fcdf(Fstat,dofr-doff,doff,'upper');
        end
        [p_fdr,p_masked]=fdr(pval{s},alpha);
        issig{e,s}=p_masked;
        
    end
    
end

%%
r=cellfun(@mean,issig);
r2=cellfun(@sum,issig);

%%
figure
hs(1)=subplot(311);
bar(1:20,100*r');
set(get(hs(1),'ylabel'),'String','% significant voxels');
set(get(hs(1),'xlabel'),'String','Subject');
hlg0=legend('Echo 1','Echo 2','Echo 3'); 

%
%figure
% find best subjects
[~,sortind]=sort(mean(r),'descend');
nSubjShow=8;
TR=2.8;
time=(0:644)*TR;
for s=1:nSubjShow
    hs(s+1)=subplot(6,2,4+s); hold on
    for e=1:3
        thisTs=allts{e,sortind(s)};
        vinds=issig{e,sortind(s)};
        muThisTs=mean(thisTs(:,vinds==1),2);
        plot(time,muThisTs);
    end
    %ylim([-1 1]);
    axis tight
    yl=ylim;
    plot([600 600],[yl(1) yl(2)],'-k');
    plot([1200 1200],[yl(1) yl(2)],'-k');
    harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
    
    %htxt=text(50,yl(2)-0.25,['S' num2str(sortind(s))]);
    htxt=text(50,yl(1)+0.15,['S' num2str(sortind(s))]);
    
    if s==1, 
        hlg=legend('Echo 1','Echo 2','Echo 3'); 
        htxt=text(725,yl(1)+0.25,['Light on']);
    end
    
end

hins = axes('position', [0.14 0.85 0.2 0.05]);
plot(hins,1:999,[zeros(1,333) ones(1,333) zeros(1,333)],'k','LineWidth',2); box off
htxt=text(360,1.37,'Design');
set(hins,'XTick',[]); set(hins,'YTick',[]);

% esthetics
set(hs(1),'Xtick',1:20);
set(hs(1),'Box','off');

set(hs(2:end),'XLim',[0 1800]);
set(hs(2:end),'XTick',[0 600 1200 1800]);
set(hs(2:nSubjShow-1),'XTickLabel',{'','','',''});
set(hs(nSubjShow),'XTickLabel',{'0:00','10:00','20:00','30:00'});
set(hs(nSubjShow+1),'XTickLabel',{'0:00','10:00','20:00','30:00'});

set(hs(2:end), 'YLimSpec', 'Tight');

%set(hs(2:end),'Ytick',[-1 0 1]);
%set(hs(3:2:end),'YTickLabel',{'','',''});

set(get(hs(2),'YLabel'),'String','BOLD (a.u.)');
set(get(hs(nSubjShow),'XLabel'),'String','Time (mm:ss)');

set(hlg0,'box','off','Orientation','horizontal');

lgPos=get(hlg,'Position');
delx=-0.07; dely=-0.07;
set(hlg,'Position',[lgPos(1)+delx lgPos(2)+dely lgPos(3) lgPos(4)]);
set(hlg,'box','off','Orientation','horizontal');

sublabel([hs(1) hs(2)],-10,-60,'FontWeight','Bold','FontSize',16);

if PRINT_FIGURE
print -dpng ../figures/revfig1
end

% % %%
s=2;
b1=allBetas{1,s}(issig{1,s});
b2=allBetas{2,s}(issig{1,s});
b3=allBetas{3,s}(issig{1,s});
figure;
hold on
plot(ones(1,numel(b1)),b1,'o');
plot(2*ones(1,numel(b2)),b2,'o');
plot(3*ones(1,numel(b3)),b3,'o');