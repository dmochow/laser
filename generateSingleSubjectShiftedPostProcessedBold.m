clear all; close all; clc

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);
CSHIFTOPTION=2;
cshifts=computeCshifts(CSHIFTOPTION);

%%
for s=1:nSubjects
    for e=1:3
        [s,e]
        ECHOSTR=num2str(e);
        pathToData=['../data/' subjStrs{s} '/NII'];
        boldFilename=['s8oBoldEcho' ECHOSTR '+tlrc.BRIK'];
        [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,boldFilename));
        shBold=circshift(bold,[cshifts(1,s) cshifts(2,s) cshifts(3,s) 0]);
        [err,ErrMessage,Info]=WriteBrikWrap(pathToData,shBold,boldInfo,['shs8oBoldEcho' ECHOSTR],'tlrc');
    end
end

%%
% origPath=pwd;
% for s=1:nSubjects
%     pathToData=['../data/' subjStrs{s} '/NII'];
%     cd(pathToData);
%     delete dsbold_e1_al+orig.BRIK
%     delete dsbold_e2_al+orig.BRIK
%     delete dsbold_e3_al+orig.BRIK
%     
%     delete shs8oBoldEcho3+tlrc.HEAD
%     delete shs8oBoldEcho3+tlrc.BRIK
%     
%     delete shs8oBoldEcho2+tlrc.HEAD
%     delete shs8oBoldEcho2+tlrc.BRIK
%     
%     delete shs8oBoldEcho1+tlrc.HEAD
%     delete shs8oBoldEcho1+tlrc.BRIK
%     
%     cd(origPath);
% end