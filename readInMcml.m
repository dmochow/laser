clear all; close all; clc;
%dataFilename='../data/mcml/moreTits_trimmed_v3.mat';
dataFilename='../data/mcml/revisedTLS_contourData_trimmed.mat';
load(dataFilename)
whos
%A=averageMap;
A=truncatedAveragedData; % the one I want
Z=(0:size(A,1)-1)*dz;
R=(0:size(A,2)-1)*dr;
Z=Z*10; % convert to mm
R=R*10;
save('../data/precomputed/RZAfull_092418.mat','R','Z','A');
%save('../data/precomputed/RZAfull.mat','R','Z','A');