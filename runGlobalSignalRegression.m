% regress out global signal
clear all; close all; clc
precomputedFilename=['../data/precomputed/allBoldsRoi ' date];
load(precomputedFilename,'allBoldsRoi','allBoldsGlobal','allOnsets','allBoldsControlRoi','TR','TEs','subjStr','nSubjects','nEchos');

for s=1:nSubjects
    for e=1:nEchos
        thisBold=squeeze(allBoldsRoi{s}(e,:,:));
        thisBoldGlobal=squeeze(allBoldsGlobal{s}(e,:)).';
        %tb=tbg*a
        %a=pinv(tbg)*tb
        %thisBoldOut=transpose( thisBold.'-thisBold.'*pinv(thisBoldGlobal)*thisBoldGlobal );
        %allBoldsRoiOut{s}(e,:,:)=thisBoldOut;
        allBoldsRoiOut{s,1}(e,:,:) = regressOut(thisBold,thisBoldGlobal).';
    end
end

%%
save(precomputedFilename,'allBoldsRoi','allBoldsGlobal','allOnsets','allBoldsControlRoi','TR','TEs','subjStr','nSubjects','nEchos','allBoldsRoiOut');


%%
allBolds=cat(3,allBoldsRoiOut{:});
muBold=mean(allBolds,3);

%%
figure
for e=1:nEchos
    subplot(nEchos,1,e);
    plot(muBold(e,:));
end
