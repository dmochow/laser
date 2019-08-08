% 06/05/18
% optimize the AR model order so that Chow analysis can be parameterized
% appropriately

clear all; close all; clc
TR=2.8; nTR=645; time=(0:nTR-1)*TR;

dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo1-02-Jun-2018.mat';
%dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo2-31-May-2018.mat';
%dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo3-31-May-2018.mat';

load(dataFilename,'alloBolds');as 
allBolds2D=cat(2,alloBolds{:});
grandBold=mean(allBolds2D,2);
y=grandBold;
nSamples=numel(y);

orders=[2:11];
nOrders=numel(orders);

nFolds=5;
mse=zeros(nOrders,nFolds);
for o=1:nOrders
    nlags=orders(o);
    X=tplitz(y,nlags); X=X(:,2:end-1);
    
    nSamplesPerFold=floor(nSamples/nFolds);
    nSamplesTrunc=nSamplesPerFold*nFolds;
    cvInds=reshape(1:nSamplesTrunc,nSamplesPerFold,nFolds);
    
    for f=1:nFolds
        [o,f]
        trainSamples=cvInds(:,setdiff(1:nFolds,f)); trainSamples=trainSamples(:);
        testSamples=cvInds(:,f); testSamples=testSamples(:);
        Xtrain=X(trainSamples,:); ytrain=y(trainSamples);
        Xtest=X(testSamples,:); ytest=y(testSamples);
        h=pinv(Xtrain)*ytrain;
        allH(1:nlags-1,o,f)=h;
        mse(o,f)=mean((ytest-Xtest*h).^2);
    end
    
end

%%
J=mean(mse,2);
[~,optOrder]=min(J);
%% draw figure
figure;
subplot(221);
plot(orders-1,J,'-ok','MarkerFaceColor','k');
xlabel('AR Model Order');
ylabel('Mean Squared Error');
subplot(222);
stem(1:max(orders)-1,allH(:,optOrder,1),'k');
xlabel('Lag');
ylabel('AR Coefficient');
xlim([0 11]);
print -dpng ../figures/arModelOrder
