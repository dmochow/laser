% 04/23/18
% toy var
clear all; close all; clc

D=3;
p=1;
T=5000;
%A=randn(D,D,p);  % AR matrix
A=[1,-0.7,0.3 ; 1,0.5,-0.2; 1,-0.2,0.1];
e=0.1*randn(T,D); 
numobs=T; % number of samples to forecast
numpaths=5;

% generate some synthetic data
Y=zeros(T,D);
Y(1:p,:)=zeros(p,D); % initial p values
for t=p+1:T
    for l=1:p
        Y(t,:)=Y(t,:)+Y(t-l,:)*A(:,:,l); 
    end
    Y(t,:)=Y(t,:)+e(t,:);
end


Yf = fitAndForecast(Y,p,numobs,numpaths);

% % create model
% myAR=varm(D,p);
% 
% % fit model
% myARest=estimate(myAR,Y);
% 
% % now run forecast
% Yf = simulate(myARest,numobs,'Y0',Y(end-p+1:end,:),'NumPaths',numpaths);

figure; 
subplot(311); plot([Y; Yf(:,:,1) ]);
subplot(312); plot([Y; Yf(:,:,2) ]);
subplot(313); plot([Y; Yf(:,:,3) ]);
