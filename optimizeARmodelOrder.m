% 09.29.18 trying only on grey matter
% 06/05/18
% 09.24.18
% need to determine: optimal order for each echo, and whether one should or
% should not include the constant term
% thus we need to construct 6 matrices of P x F x V


% no constant minJ = 3.2988e-04 for 100 voxels
% with constant minJ = 3.3072e-04 for 100 voxels

% no constant 5.7148e-04 random 100
% with constant 5.7277e-04 random 100


% optimize the AR model order so that Chow analysis can be parameterized
% appropriatelyf
clear all; close all; clc
COMPUTE_ALL=1; % actually run the exhaustive search (or just show the result)
nEchos=3; TR=2.8; nTR=645; time=(0:nTR-1)*TR;
pathToData='../data/SAVG/NII/';
brainMaskFilename='resampled_brain_mask+tlrc';
orders=[2:11]; % p=orders-1
nOrders=numel(orders);
nFolds=5;
CONSTANT=1;
nSamplesPerFold=floor(nTR/nFolds);
nSamplesTrunc=nSamplesPerFold*nFolds;
cvInds=reshape(1:nSamplesTrunc,nSamplesPerFold,nFolds);

[~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
brainMask=logical(brainMask);


% %% resample segmentation onto bold
%    !export PATH=$PATH:/Users/jacekdmochowski/abin
%     PATH = getenv('PATH');
%     setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
%     origPath=pwd;
%     cd(pathToData);
%     str=['!3dresample -master ' brainMaskFilename ' -prefix resampled_Classes -input Classes+tlrc'];
%     eval(str);
%     cd(origPath);
% %%    
    
[~,labels,~,~]=BrikLoad('../data/SAVG/NII/resampled_Classes+tlrc.BRIK');
isGM=labels==2;

gmMask=brainMask & isGM;
nVoxels=sum(gmMask(:));

% get all data into memory
smuBold2D=zeros(nTR,nVoxels,nEchos);
for e=1:nEchos
    ECHOSTR=num2str(e);
    dataFilename=['s8muwBoldEcho' ECHOSTR '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc'];
    [~, smuBold, info, ~] = BrikLoad (fullfile(pathToData,dataFilename));
    %smuBold2D(:,:,e)=vol2ts(smuBold,brainMask);
    smuBold2D(:,:,e)=vol2ts(smuBold,gmMask);
end


% quick runs (uncomment next 3 lines)
%perms=randperm(nVoxels);
%nVoxels=800;
%smuBold2D=smuBold2D(:,perms(1:nVoxels),:);

%%
if COMPUTE_ALL
    
    mse=zeros(nOrders,nFolds,nVoxels,nEchos);
    
    %%
    for e=1:nEchos
        for v=1:nVoxels
            [e,v]
            y=smuBold2D(:,v,e); % the data for this voxel and echo
            y=y/sum(y.^2);
            for o=1:nOrders
                nlags=orders(o);
                X=tplitz(y,nlags);
                
                if CONSTANT
                    X=X(:,2:end);
                else
                    X=X(:,2:end-1);
                end
                
                for f=1:nFolds
                    trainSamples=cvInds(:,setdiff(1:nFolds,f)); trainSamples=trainSamples(:);
                    testSamples=cvInds(:,f); testSamples=testSamples(:);
                    Xtrain=X(trainSamples,:); ytrain=y(trainSamples);
                    Xtest=X(testSamples,:); ytest=y(testSamples);
                    h=pinv(Xtrain)*ytrain;
                    mse(o,f,v,e)=mean((ytest-Xtest*h).^2);
                end
                
            end
        end
    end
    %%
    J=nanmean(nanmean(nanmean(mse,4),3),2);
    [minJ,optOrder]=min(J);
    
    
    
    
    %% draw figure
    figure;
    subplot(221);
    plot(orders-1,J,'-ok','MarkerFaceColor','k');
    xlabel('AR Model Order');
    ylabel('Mean Squared Error');
    print('-dpng',['../figures/arModelOrderGMonlyConstant' num2str(CONSTANT) '.png']);
    save(['../data/precomputed/arModelOptimizationGMonlyConstant' num2str(CONSTANT)],'orders','mse','J','minJ','optOrder');
end

%% combine both sets of results (w and wo constant term)
load('../data/precomputed/arModelOptimizationGMonlyConstant1','orders','mse','J','minJ','optOrder');
J1=J; mse1=mse;
load('../data/precomputed/arModelOptimizationGMonlyConstant0','orders','mse','J','minJ','optOrder');
J0=J; mse0=mse;
[orderKnee1, orderKneeIndx1] = knee_pt(J1,orders); % p=5;
[orderKnee0, orderKneeIndx0] = knee_pt(J0,orders); % 
nJ0=J0/max(J0);
nJ1=J1/max(J1);
figure;
hs(1)=subplot(221);
plot(orders-1,nJ0,'k'); hold on
hk=plot(orders(orderKneeIndx0)-1,nJ0(orderKneeIndx0),'ob','MarkerFaceColor','b');
title('Without constant term','FontWeight','normal');
ylabel('Normalized Error');
xlabel('AR order');
hlg=legend([hk],'Knee Point');
set(hlg,'box','off');

hs(2)=subplot(222);
plot(orders-1,nJ1,'k'); hold on
plot(orders(orderKneeIndx1)-1,nJ1(orderKneeIndx1),'ob','MarkerFaceColor','b');
title('With constant term','FontWeight','normal');
xlabel('AR order');
hlg=legend([hk],'Knee Point');
set(hlg,'box','off');
sublabel(hs,-20,-20,'FontWeight','bold','FontSize',16);
print('-dpng',['../figures/arModelOrder.png']);
crop('../figures/arModelOrder.png',0);
