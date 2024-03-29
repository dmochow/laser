clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

% compute the group mean
str=['!3dmean -prefix ../data/SAVG/groupMean_e1 ' ...
    '../data/S02/NII/smnudsbold_e1_tlrc_al+tlrc ' ... 
    '../data/S03/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S04/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S05/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S06/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S07/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S09/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S10/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S11/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S12/NII/smnudsbold_e1_tlrc_al+tlrc ' ...
    '../data/S13/NII/smnudsbold_e1_tlrc_al+tlrc '];
eval(str);

% echo 2
str=['!3dmean -prefix ../data/SAVG/groupMean_e2 ' ...
    '../data/S02/NII/smnudsbold_e2_tlrc_al+tlrc ' ... 
    '../data/S03/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S04/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S05/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S06/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S07/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S09/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S10/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S11/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S12/NII/smnudsbold_e2_tlrc_al+tlrc ' ...
    '../data/S13/NII/smnudsbold_e2_tlrc_al+tlrc '];
eval(str);

% echo 3
str=['!3dmean -prefix ../data/SAVG/groupMean_e3 ' ...
    '../data/S02/NII/smnudsbold_e3_tlrc_al+tlrc ' ... 
    '../data/S03/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S04/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S05/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S06/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S07/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S09/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S10/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S11/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S12/NII/smnudsbold_e3_tlrc_al+tlrc ' ...
    '../data/S13/NII/smnudsbold_e3_tlrc_al+tlrc '];
eval(str);