% smooth individually pre- and post-processed bolds in preparation for chow
% testing on individual level
clear all; close all; clc

!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

FWHMSTR='8';

nEchos=3;

%['oBoldEcho' ECHOSTR],'tlrc'

for s=1:nSubjects
    for e=1:nEchos
    
    inputFilename=['oBoldEcho' num2str(e) '+tlrc'];
    prefix=cat(2,'s8',inputFilename);    
    
    origPath=pwd;
    cd(['../data/' subjStrs{s} '/NII/']);
    
    switch e
        case 1
            delete s8oBoldEcho1+tlrc.BRIK
            delete s8oBoldEcho1+tlrc.HEAD
        case 2
            delete s8oBoldEcho2+tlrc.BRIK
            delete s8oBoldEcho2+tlrc.HEAD
        case 3
            delete s8oBoldEcho3+tlrc.BRIK
            delete s8oBoldEcho3+tlrc.HEAD
    end
    
    %muwPrefix=['s' tFWHMSTR 'muwBoldEcho' optionStr];
    str=['!3dBlurInMask -prefix ' prefix ' -input ' inputFilename  ' -FWHM ' FWHMSTR ' -mask ' inputFilename];
    eval(str);
    cd(origPath);
    
    end
end