function [coeffs,arMagSpec,pvals,tstats] = computeArCoeffs(groupBoldTs)
% for a given (group) bold, compute the AR coefficients and spectra
% groupBoldTs is size nTRs by nVoxels
NLAGS=5; %p=NLAGS-1 
LASERONSETTR=215;
LASEROFFSETTR=430;
nfft=32;

%%
nVoxels=size(groupBoldTs,2);
coeffs=zeros(NLAGS-1,nVoxels,3); 
pvals=zeros(nVoxels,2);
tstats=zeros(nVoxels,2);

for v=1:nVoxels
    v
    y=groupBoldTs(:,v);
    X=tplitz(y,NLAGS);
    X=X(:,2:end-1); % don't include constant term
    
    %% pre versus stim
    [~,p1,tstat1] = chowtest(X(1:LASEROFFSETTR,:),y(1:LASEROFFSETTR),LASERONSETTR);
    [~,p2,tstat2] = chowtest(cat(1,X(1:LASERONSETTR,:),X(LASEROFFSETTR+1:end,:)),cat(1,y(1:LASERONSETTR),y(LASEROFFSETTR+1:end)),LASERONSETTR);
    
    X1=X(1:LASERONSETTR,:); y1=y(1:LASERONSETTR);
    X2=X(LASERONSETTR+1:LASEROFFSETTR,:); y2=y(LASERONSETTR+1:LASEROFFSETTR);
    X3=X(LASEROFFSETTR+1:end,:); y3=y(LASEROFFSETTR+1:end);
    
    coeffs(:,v,1)=X1\y1;
    coeffs(:,v,2)=X2\y2;
    coeffs(:,v,3)=X3\y3;
        
    pvals(v,1)=p1; pvals(v,2)=p2;
    tstats(v,1)=tstat1;  tstats(v,2)=tstat2; 
    
    
    
end
arMagSpec = getARspectrum(coeffs,nfft);

