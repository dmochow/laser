clear all; %close all; clc

% 
%slice=5;
%vox=[113,136; 120,100; 128,61];
subjStr='tS01';
% pathToData=['../data/thermo/' subjStr '/EPI_THERM_WIP1118_3D_PHASE_CHANGE_0022/'];
% pathToRefData=['../data/thermo/' subjStr '/EPI_THERM_WIP1118_3D_PHASE_CHANGE_0023/'];
% pathToMagData=['../data/thermo/' subjStr '/EPI_THERM_WIP1118_3D_MAG_0020/'];

pathToData=['../data/thermo/' subjStr '/'];
pathToRefData=['../data/thermo/' subjStr '/'];
pathToMagData=['../data/thermo/' subjStr '/'];


% get dicom info
%files=dir(fullfile(pathToData,'*.ima'));
%info = dicominfo(fullfile(pathToData,files(1).name));
delrad=8*pi/(2^16-1);
%slope=info.RescaleSlope;
slope=delrad;
%intercept=info.RescaleIntercept;
intercept=-12.56;

%nii=load_untouch_nii(fullfile(pathToData,'data.nii'));
nii=load_untouch_nii(fullfile(pathToData,'ph.nii'));
x=nii.img;
x=double(intercept+double(x)*slope);

nii=load_untouch_nii(fullfile(pathToRefData,'phb.nii'));
xr=nii.img;
xr=double(intercept+double(xr)*slope);
%xr=double(xr);

nii=load_untouch_nii(fullfile(pathToMagData,'mag.nii'));
anat=nii.img;
anat=double(anat);

%%
xr4=repmat(xr,[1 1 1 size(x,4)]);

%%
% convert to temperature change
a=-0.01/1e6; % "a" Yuan et al. 2012
gamma=42.58e6; % gyromagnetic ratio of hydrogen in water
Bo=3; % 3 Tesla
TE=17e-3; % 17 ms?

%y=(x-xr4)/(a*gamma*Bo*TE);
y=(x)/(a*gamma*Bo*TE);

%%
% display slices
figure;
sl=[2:2:32]; 
TR=round(20+3*(150/20.5));
%img=squeeze(y(:,:,sl,TR));
for s=1:numel(sl)
    hs(1)=subplot(3,6,s);
    imagesc(squeeze(y(:,:,sl(s),TR)));
end
colormap jet

%%
% display anat slices
figure;
sl=[2:2:32]; 
TR=round(20+3*(150/20.5));
%img=squeeze(y(:,:,sl,TR));
for s=1:numel(sl)
    hs(1)=subplot(3,6,s);
    imagesc(squeeze(anat(:,:,sl(s),TR)));
end
colormap jet

%%
% display time courses
ts1=squeeze(y(39,53,16,:));
ts2=squeeze(y(38,52,16,:));
ts3=squeeze(y(26,51,16,:));
ts4=squeeze(y(27,50,16,:));
figure; hold on; plot(ts1); plot(ts2); plot(ts3); plot(ts4);
legend('L','L','R','R');