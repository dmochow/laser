function oBold = smoothBold(pathToData,filename,fhwm,orig)
%a wrapper that lets you smooth afni briks without you yourself having to
%go to the command line
%
% TODO: option to get BRIK file written
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

if nargin<4, orig=0; end % default to tlrc

if orig==1
    oPrefix='tmp+orig';
else
    oPrefix='tmp+tlrc';
end
    

origPath=pwd;
cd(pathToData);

stro=['delete ' oPrefix '.BRIK']; eval(stro);
stro=['delete ' oPrefix '.HEAD']; eval(stro);

str=['!3dBlurInMask -prefix ' oPrefix ' -input ' filename ' -FWHM ' num2str(fhwm) ' -mask ' filename];
eval(str);

[~, oBold, ~, ~] = BrikLoad (oPrefix); % we're already in pathToData

cd(origPath);
        
end

