clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%% not working due to stdin problem
str='!3dDeconvolve -nodata 300 1 -polort -1 -num_stimts 2 -stim_times 1 ''1D: 10 150'' ''MION(70)'' -stim_times 2 ''1D: 10 150'' ''BLOCK(70,1)'' -x1D stdout: | 1dplot -stdin -one -thick  ';
eval(str);


%%
3dDeconvolve -nodata 300 1  -num_stimts 2   \
                            -stim_times 1 '1D: 10 150' 'MION(70)'    \
                            -stim_times 2 '1D: 10 150' 'BLOCK(70,1)' \
                            -x1D stdout: | 1dplot -stdin -one -thick   