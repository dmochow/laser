clear all; close all; clc

% base parameters
t=linspace(0,0.1,1000);
T2s=0.035;
So=1;
TEs=[12.8,34.3,55.6]/1000;

% simulated bold
S=So*exp(-t/T2s);
for e=1:numel(TEs), [~,TEsmp(e)]=min(abs(t-TEs(e))); end

% simulated change in T2 (increase so longer time constant)
delT2s=0.1*T2s;
S2=So*exp(-t/(T2s+delT2s));

% increase in So
delSo=0.1*So;
S3=(So+delSo)*exp(-t/T2s);

% increase in So + decrease in T2s
S4=(So+delSo)*exp(-t/(T2s+delT2s));

% 
figure;
hs(1)=subplot(221);
plot(t*1000,S,'k','LineWidth',2);
hold on
plot(t*1000,S2,'Color',[1 0 0]*0.75,'LineWidth',2);
hlg=legend('baseline','T2^\ast increase');
xlabel('Echo time (ms)','FontSize',14);
ylabel('Signal intensity','FontSize',14);
set(hlg,'box','off');
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)+0.02 lgPos(2)+0.025 lgPos(3) lgPos(4)]);
box off

hs(2)=subplot(222); hold on
plot(t*1000,(S2-S)./S*100,'Color',[1 0 0]*0.75,'LineWidth',2);
ylabel('Percent change','FontSize',14);
xlabel('Echo time (ms)','FontSize',14);

print -dpng ../figures/asrcT2_2of2
crop('../figures/asrcT2_2of2.png');

figure;
hs(1)=subplot(221);
plot(t*1000,S,'k','LineWidth',2);
hold on
%plot(t*1000,S2,'Color',[1 0 0]*0.75,'LineWidth',2);
hlg=legend('baseline');
xlabel('Echo time (ms)','FontSize',14);
ylabel('Signal intensity','FontSize',14);
set(hlg,'box','off');
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)+0.02 lgPos(2)+0.025 lgPos(3) lgPos(4)]);
box off

hs(2)=subplot(222); hold on
%plot(t*1000,(S2-S)./S*100,'Color',[1 0 0]*0.75,'LineWidth',2);
ylabel('Percent change','FontSize',14);
xlabel('Echo time (ms)','FontSize',14);

print -dpng ../figures/asrcT2_1of2
crop('../figures/asrcT2_1of2.png');



figure;
hs(1)=subplot(221);
plot(t*1000,S,'k','LineWidth',2);
hold on
plot(t*1000,S3,'Color',[0 1 0]*0.75,'LineWidth',2);
% plot(t*1000,S4,'Color',[0 1 0]*0.75,'LineWidth',2);
hlg=legend('baseline','S_o increase');
xlabel('Echo time (ms)','FontSize',14);
ylabel('Signal intensity','FontSize',14);
set(hlg,'box','off');
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)+0.02 lgPos(2)+0.025 lgPos(3) lgPos(4)]);
box off

hs(2)=subplot(222); hold on
plot(t*1000,(S3-S)./S*100,'Color',[0 1 0]*0.75,'LineWidth',2);
ylabel('Percent change','FontSize',14);
xlabel('Echo time (ms)','FontSize',14);
ylim([0 30])

print -dpng ../figures/asrcSo_2of2
crop('../figures/asrcSo_2of2.png');

figure;
hs(1)=subplot(221);
plot(t*1000,S,'k','LineWidth',2);
hold on
%plot(t*1000,S3,'Color',[0 1 0]*0.75,'LineWidth',2);
hlg=legend('baseline');
xlabel('Echo time (ms)','FontSize',14);
ylabel('Signal intensity','FontSize',14);
set(hlg,'box','off');
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)+0.02 lgPos(2)+0.025 lgPos(3) lgPos(4)]);
box off
ylim([0 1.2]);

hs(2)=subplot(222); hold on
%plot(t*1000,(S3-S)./S*100,'Color',[0 1 0]*0.75,'LineWidth',2);
ylabel('Percent change','FontSize',14);
xlabel('Echo time (ms)','FontSize',14);
ylim([0 30])

print -dpng ../figures/asrcSo_1of2
crop('../figures/asrcSo_1of2.png');


%%
figure;
hs(1)=subplot(221);
plot(t*1000,S,'k','LineWidth',2);
hold on
plot(t(300)*1000,S(300),'ko','LineStyle','none','MarkerFaceColor','k');
xlabel('Echo time (ms)','FontSize',14);
ylabel('Signal intensity','FontSize',14);
box off

hs(2)=subplot(222); hold on
plot(t*1000,S,'k','LineWidth',2);
hold on
plot(t(100)*1000,S(100),'ko','LineStyle','none','MarkerFaceColor','k');
plot(t(300)*1000,S(300),'ko','LineStyle','none','MarkerFaceColor','k');
plot(t(500)*1000,S(500),'ko','LineStyle','none','MarkerFaceColor','k');
xlabel('Echo time (ms)','FontSize',14);
box off

print -dpng ../figures/asrc_mebold
crop('../figures/asrc_mebold.png');


% 
% sublabel(hs,0,-40,'FontWeight','Bold','FontSize',16);
% print -dpng ../figures/simulatedEffect
% crop('../figures/simulatedEffect.png');
% 
% % subplot(222);
% % 
% % hold on
% % plot(t(TEsmp),S2(TEsmp),'or');
% 

