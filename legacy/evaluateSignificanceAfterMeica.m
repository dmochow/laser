clear all; close all; clc
% take prantik's output and determine if any voxels show effect of laser
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));
subjIndx=4;
alpha=0.05; % for output nii
nPerm=500; % for significance testing
%%
switch subjIndx
    case 2
        maskFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/meica.bold_e123/eBvrmask.nii'];
        dataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn.nii'];
        biopacFilename=['../data/S0' num2str(subjIndx) '/BIOPAC/greg.mat'];
        outFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn_pvals.nii'];
        saveDataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/matData'];
        TR=2.8;
        fs=1/2.8;
        DUR_frames=645; % duration of BOLD recording in frames
    case 3
        maskFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/meica.bold_e123/eBvrmask.nii'];
        dataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn.nii'];
        biopacFilename=['../data/S0' num2str(subjIndx) '/BIOPAC/subj0003_051917.mat'];
        outFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn_pvals.nii'];
        saveDataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/matData'];
        TR=2.8;
        fs=1/2.8;
        DUR_frames=645; % duration of BOLD recording in frames
     case 4
        maskFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/meica.bold_e123/eBvrmask.nii'];
        dataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn.nii'];
        biopacFilename=['../data/S0' num2str(subjIndx) '/BIOPAC/subj0004_053117.mat'];
        outFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn_pvals.nii'];
        saveDataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/matData'];
        TR=2.8;
        fs=1/2.8;
        DUR_frames=645; % duration of BOLD recording in frames  
     case 5
        maskFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/meica.bold_e123/eBvrmask.nii'];
        dataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn.nii'];
        biopacFilename=['../data/S0' num2str(subjIndx) '/BIOPAC/subj0005_061217.mat'];
        outFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/bold_e123_medn_pvals.nii'];
        saveDataFilename=['../data/S0' num2str(subjIndx) '/NII/mebold2go/matData'];
        TR=2.8;
        fs=1/2.8;
        DUR_frames=645; % duration of BOLD recording in frames    
end

%% figure out on/off times of laser from biopac
load(biopacFilename);
trigger=data(:,1); % discard second channel
diffTrigger=diff(trigger);
[~,timeOn_ms]=max(diffTrigger);
[~,timeOff_ms]=min(diffTrigger);
TR_ms=TR*1000;
timeOn_frames=round(timeOn_ms/TR_ms);
timeOff_frames=round(timeOff_ms/TR_ms);
tOn=timeOn_frames:timeOff_frames-1;
tOff_1=1:timeOn_frames-1;
tOff_2=timeOff_frames:DUR_frames;
tOff=[tOff_1 tOff_2];
% show the result to check correctness
figure
subplot(211);
plot( (1:numel(diffTrigger)) / (1000*TR ) ,  diffTrigger);
yl=ylim;xl=xlim;
subplot(212);
plot(timeOn_frames*[1 1],linspace(yl(1),yl(2),2),'--r');hold on
plot(timeOff_frames*[1 1],linspace(yl(1),yl(2),2),'--b');xlim(xl);
legend('estimated on','estimated off'); legend boxoff;
%%
maskNii=load_untouch_nii(maskFilename);
mask=logical(maskNii.img);

%%
dataNii=load_untouch_nii(dataFilename);
data=double(dataNii.img);
data_4D=data;
data=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
brainData=data(mask(:),:);

%%
% 38,68,29
roi_x=33:43;
roi_y=63:73;
roi_z=25:35;

roiData=squeeze(data_4D(roi_x,roi_y,roi_z,:));
roiData2=permute(roiData,[4 1 2 3]);
roiData3=roiData2(:,:);
muRoiData=mean(roiData3,2);
figure; plot(muRoiData,'k');
%%
figure; plot(smoothts(muRoiData,'g',500));

%%
figure; imagesc(data_4D(:,:,30)); axis xy; colormap bone

%%
if 1
    
    
    
    %%
    % for each voxel in the brain, compute fit of time series to the design
    % function (flat-ramp-flat)
    %%
    design=zeros(numel(tOff_1)+numel(tOn)+numel(tOff_2),1);
    design(tOn)=linspace(0,1,numel(tOn));
    design(tOn)=1;
    %design(tOff_2)=linspace(1,0,numel(tOff_2));
    nBrainVoxels=numel(find(mask(:)));
    r2=zeros(nBrainVoxels,1);
    r2_mock=zeros(nBrainVoxels,nPerm);
    
    
    
    Y=brainData;
    X=repmat(design(:)',nBrainVoxels,1);
    XY=sum(X.*Y,2);
    XX=sum(X.*X,2);
    BETA=XY./XX;
    BETA=repmat(BETA,1,DUR_frames);
    sse=sum((Y-BETA.*X).^2,2);
    
    % permutation test
    sse_mock=zeros(size(sse,1),nPerm);
    for p=1:nPerm
        p
        Y_mock=(surrogateResponseGenerator(Y.')).';
        XY_mock=sum(X.*Y_mock,2);
        BETA_mock=XY_mock./XX;
        BETA_mock=repmat(BETA_mock,1,DUR_frames);
        sse_mock(:,p)=sum((Y_mock-BETA_mock.*X).^2,2);
    end
    
    %%
    SSE=repmat(sse,1,nPerm);
    pvals=mean(SSE>sse_mock,2);
    
    % correct p-values for multiple comparisons
    [p_fdr, p_masked] = fdr( pvals, alpha);
    
    % %
    tmp=zeros(size(mask(:)));
    %tmp(mask(:))=p_masked;
    tmp(mask(:))=1-pvals;
    pvalVolume=reshape(tmp,[size(mask,1),size(mask,2),size(mask,3)]);
    outNii=maskNii;
    outNii.img=pvalVolume;
    save_untouch_nii(outNii,outFilename);
    %
    % save(saveDataFilename,'r2','r2_mock','pvals','p_fdr','p_masked','pvalVolume');
    
end


%
% %%
% % basic figure
% clrs=[0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250];










