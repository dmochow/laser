% contour plot of absorption data
clear all; close all; clc

load ../data/mcml/forJacek_trimmed

xAxis = linspace(0,axialCutOff, size(thinnedTruncatedAverageData,1));
yAxis = linspace(0,radialCutOff, size(thinnedTruncatedAverageData,1));

hLines = [0.6, 1.1, 1.2];

[X,Y] = meshgrid(xAxis, yAxis);
figure;
hold on; box on;
contour(X,Y,thinnedTruncatedAverageData, [0.1 0.5 1 10 100], 'ShowText', 'on')
xlabel('Radial distance (cm)');
ylabel('Axial distance (cm)');
legend('J/cm^3');

contour(-X,Y,thinnedTruncatedAverageData, [0.1 0.5 1 10 100], 'ShowText', 'on')
xlabel('Radial distance (cm)');
ylabel('Axial distance (cm)');
legend('Energy Absorbed (J/cm^3)');
line([-3, 3],[0.6, 0.6], 'Color', 'red', 'LineStyle', '--');
line([-3, 3],[1.1, 1.1], 'Color', 'red', 'LineStyle', '--');
line([-3, 3],[1.2, 1.2], 'Color', 'red', 'LineStyle', '--');
xlim([-2,2]);
ylim([0, 2.5]);
set(gca, 'fontsize', 20)