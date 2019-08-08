clear all; close all; clc

load ../data/mcml/Arz

x=[r(1) r(end)];
y=[-z(end) z(end)];

Adraw=cat(2,Abs(end:-1:1,:).',Abs.');

figure;
imagesc(y,x,Adraw);
%surf(Abs);
ylabel('depth (cm)')
xlabel('radius (cm)');
xlim([-2 2]); % 2 cm radial
ylim([0 3]); % 3 cm depth
title('Absorption (J/cm^2)','FontWeight','normal','FontSize',16);
hcb=colorbar;
colormap(jet);
print -dpng ../figures/mcml

%% 
% determine the ROI here
thresh=0.99;
figure;
subplot(221)
cs=cumsum(sum(Abs,1));
ncs=cs/max(cs);
plot(z, ncs );
zo=z(find(ncs>thresh,1));

subplot(222)
cs=cumsum(sum(Abs,2));
ncs=cs/max(cs);
plot(r, ncs );
ro=r(find(ncs>thresh,1));

[zo,ro]