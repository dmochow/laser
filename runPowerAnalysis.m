% 06/01/18
% run power analysis for TRISH
%
clear all; close all; clc
TR=2.8; nTR=645;
time=(0:nTR-1)*TR;
dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-TSHIFT1-Echo1-05-Jun-2018'
load(dataFilename,'alloBolds');
allBolds2D=cat(2,alloBolds{:});
nReps=320; % from protocol
%%
% power analysis here

tm1=780;
tm0=551.6;

[~,t1]=min(abs(time-tm1)); % treatment
[~,t0]=min(abs(time-tm0)); % control

x1=allBolds2D(t1,:); x1=x1(:);
x0=allBolds2D(t0,:); x0=x0(:);

mu1=mean(x1);
mu0=mean(x0);
std1=std(x1);
std0=std(x0);
stdp=0.5*(std1+std0);

dp=(mu1-mu0)/stdp;

%
testtype='t2';
p0=[mu0 std0];
p1=mu1;
%nout = sampsizepwr(testtype,p0,p1,pwr)
ss=[2:1:40];
nss=numel(ss);
for s=1:nss
    pwr(s) = sampsizepwr(testtype,p0,p1,[],ss(s)*nReps);
end

figure;

%
tLaser=[600 1200];
[sems,mus]=nansem(allBolds2D,2);
%figure;
hs(1)=subplot(221); hold on
time=(0:nTR-1)*TR;
shadedErrorBar(time,mus,sems,'k');
yl=ylim;
plot([tm0 tm0],[yl(1) yl(2)],'--k');
plot([tm1 tm1],[yl(1) yl(2)],'--k');
plot(time(t1),mus(t1),'or','MarkerFaceColor','r');
plot(time(t0),mus(t0),'ob','MarkerFaceColor','b');
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.75 0.75 0.75]);
set(harea,'FaceAlpha',0.25);
set(harea,'EdgeColor','none');
xlim([time(1) time(end)]);
xlabel('Time (s)');
ylabel('BOLD (a.u.)');
hlg=legend(harea,'Laser On');

hs(2)=subplot(222);
plot(ss,pwr,'k'); hold on
plot(ss(find(pwr>0.945,1)),pwr(find(pwr>0.945,1)),'ok');
%plot(ss(ss==26),pwr(ss==26),'ob');
xlabel('Sample Size N');
ylabel('Power');

sublabel(hs,0,-40,'FontSize',16,'FontWeight','Bold');
print -dpng ../figures/powerAnalysis
crop('../figures/powerAnalysis.png');

%%%%
%%
% mu1=1.14;
% mu0=1.37;
% std0=(1.48-1.25)/2;
% std1=(1.48-1.25)/2;
mu1=28.9;
mu0=28.1;
std0=0.375;
std1=0.3;
stdp=0.5*(std1+std0);
%
testtype='t2';
p0=[mu0 std0];
p1=mu1;
ss2=[2:1:20];
nss2=numel(ss2);
for s=1:nss2
    pwr2(s) = sampsizepwr(testtype,p0,p1,[],ss2(s));
end

figure;
hs(1)=subplot(221); hold on
hp=plot([1 2],[mu0;mu1],'ok','LineStyle','none');
hbar=errorbar([1 2],[mu0;mu1],[std0;std1],'LineStyle','none');
set(hbar,'Color','k');
xlim([0.5 2.5]);
set(hs(1),'XTick',[1 2]);
set(hs(1),'XTickLabel',{'Control','Treatment'});
ylabel('Correct Trials');

hs(2)=subplot(222);
plot(ss2,pwr2,'k'); hold on
plot(ss2(find(pwr2>0.945,1)),pwr2(find(pwr2>0.945,1)),'ok');
%plot(ss(ss==26),pwr(ss==26),'ob');
xlabel('Sample Size N');
ylabel('Power');

sublabel(hs,0,-40,'FontSize',16,'FontWeight','Bold');
print -dpng ../figures/powerAnalysisAim2
crop('../figures/powerAnalysisAim2.png');