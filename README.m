% README for Duke
%
% NB: afni must be installed on your machine
%
% required folder structure (sorry):

% assume that you are working in:
% WORKING_PATH/code
%
% then you must also have the following folders:
% WORKING_PATH/data/thermo/scanner_data/tSXX/
    % this folder must have
    % anat.nii
    % mag.nii
    % ph.nii
    % drift.txt
% and the preprocessed output will go in: 
% WORKING_PATH/data/thermo/output/tSXX

% prefix "t" is for thermo
% "S" is for subject
% "XX" is the subject number (e.g. 01)

% current MR-Thermo sequence:
% (1) preprocessThermoSubject.m
% (2) analyzeThermoSubject.m