% 06/05/18
% this script takes an unsmoothed post-processed bold of the form
% muwBoldEchoX and smoothes it by a bunch of different values

clear all; close all; clc
ECHOSTR='1'; % choose this to what you wany
fwhms=[8 10 12 15 20 25];
nFwhms=numel(fwhms);
for f=1:nFwhms
    FWHMSTR=num2str(fwhms(f));
    origPath=pwd;
    cd('../data/SAVG/NII/');
    delStr1=['delete s' FWHMSTR 'muwBoldEcho' ECHOSTR '+tlrc.HEAD'];
    delStr2=['delete s' FWHMSTR 'muwBoldEcho' ECHOSTR '+tlrc.BRIK'];
    eval(delStr1);
    eval(delStr2);
    str=['!3dBlurInMask -prefix s' FWHMSTR 'muwBoldEcho' ECHOSTR ' -input muwBoldEcho' ECHOSTR '+tlrc -FWHM ' FWHMSTR ' -mask muwBoldEcho' ECHOSTR '+tlrc'];
    eval(str);
    cd(origPath);
end