function Yf = fitAndForecast(Y,p,nObs,nReps)
%
% Y: data to fit model with (must be time x dimension)
% p: model order
% nObs: number of time samples to forecast
% nReps: number of forecasts

[T,D]=size(Y); % 

% create model
myAR=varm(D,p);

% fit model
myARest=estimate(myAR,Y);

if ~isstable(myARest), error('AR model is not stable.  Ejecting'); end

% now run forecast
% Yf = simulate(myARest,nObs,'Y0',Y(end-p+1:end,:),'NumPaths',nReps);
% pick a random starting point
tStart=round(T*rand)+1;
Yf = simulate(myARest,nObs,'Y0',Y(tStart-p+1:tStart,:),'NumPaths',nReps);