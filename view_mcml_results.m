clear all; close all; clc

load ../data/mcml/Arz

x=[r(1) r(end)];
y=[z(1) z(end)];

figure;
imagesc(x,y,Abs);
%surf(Abs);
xlabel('depth (cm)')
ylabel('radius (cm)');
ylim([0 2]); % 2 cm radial
xlim([0 3]); % 3 cm depth
title('Absorption (J/cm^2)','FontWeight','normal','FontSize',16);
hcb=colorbar;

print -dpng ../figures/mcml

%% 
% determine the ROI here
figure;
subplot(221)
cs=cumsum(sum(Abs,1));
ncs=cs/max(cs);
plot(z, ncs );
zo=z(find(ncs>0.99,1));

subplot(222)
cs=cumsum(sum(Abs,2));
ncs=cs/max(cs);
plot(r, ncs );
ro=r(find(ncs>0.99,1));