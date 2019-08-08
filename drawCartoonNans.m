clear all; close all; clc
% draw cartoon BOLD increase from PBM

x=1:200;
y=[randn(1,100) randn(1,100)+3];
figure;
plot(x,y,'k','LineWidth',2);
axis off
print -dpng ../figures/BOLD_cartoon.png
crop('../figures/BOLD_cartoon.png');