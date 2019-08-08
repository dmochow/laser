clear all; close all; clc
load ../data/mcml/Arz
x=[r(1) r(end)];
y=[z(1) z(end)];
figure;
subplot(221);
imagesc(Abs);
subplot(222);
imagesc(y,x,Abs);xlabel('depth'); ylabel('radius (lateral)');
subplot(223);
contour(z,r,Abs,[0.001 0.001]);

