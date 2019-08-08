function laserOrigin = fineTuneOrigin(anatFilename,laserOrigin,ulf,urf,lrf)
% project the laser onto the scalp so that we are sitting right on it
if nargin<5, error('Need 5 arguments'); end

% internal parameters
SCALP_INT_THRESH=150; % for finding edge of scalp
MAXITER=1000;
STEPSIZE=0.1; % 100 microns

initOrigin=laserOrigin; % store in case that we don't hit threshold (error)
[err, anat, Info, ErrMessage] = BrikLoad (anatFilename);

% define front face
crossprod=cross((ulf-urf),(ulf-lrf));
crossprod=crossprod/norm(crossprod);
if crossprod(2)>0 % if pointing anterior
    crossprod=-crossprod;
end
xo=laserOrigin;
imageInt=zeros(MAXITER,1);
for i=1:MAXITER
    xo=xo+crossprod*STEPSIZE;
    rxo=round(xo);
    try
        imageInt(i)=anat(rxo(1),rxo(2),rxo(3));
    catch
        imageInt(i)=imageInt(i-1);
    end
end

dImageInt=[0;diff(imageInt)];

findind=find(dImageInt>SCALP_INT_THRESH,1);
laserOrigin=laserOrigin+findind*STEPSIZE*crossprod;

% 
if ~numel(laserOrigin)
    laserOrigin=[];
    warning('Laser origin not found!');
    % the code below was not working, so doing this another way
    % by reflecting the _fine-tuned_ true laser origin
%     figure;
%     subplot(211);
%     plot(1:MAXITER,imageInt);
%     x=ginput(1);
%     maxind=x(1);
%     subplot(212);
%     plot(1:MAXITER,dImageInt);
%     laserOrigin=initOrigin+maxind*STEPSIZE*crossprod;
end

