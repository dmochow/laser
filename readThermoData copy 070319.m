clear all; close all; clc

pathToData='../data/thermo/EPI_THERM_WIP1118_3D_PHASE_CHANGE_0016/';
%pathToData='../data/thermo/EPI_THERM_WIP1118_3D_PHASE_CHANGE_0021/';
nii=load_untouch_nii(fullfile(pathToData,'data.nii'));
x=nii.img;
x=double(x);

pathToRefData='../data/thermo/EPI_THERM_WIP1118_3D_PHASE_CHANGE_0017/';
%pathToRefData='../data/thermo/EPI_THERM_WIP1118_3D_PHASE_CHANGE_0022/';
nii=load_untouch_nii(fullfile(pathToRefData,'data.nii'));
xr=nii.img;
xr=double(xr);

%%
xr4=repmat(xr,[1 1 1 size(x,4)]);

%%
% convert to temperature change
a=-0.01e6; % "a" Yuan et al. 2012
gamma=42.58; % gyromagnetic ratio of hydrogen in water
Bo=3; % 3 Tesla
TE=0.03; % 30 ms?

y=(x-xr4)/(a*gamma*Bo*TE);

yts=vol2ts(y);
xts=vol2ts(x);

%%
% figure;
% subplot(221);
% plot(mean(yts,2));
% subplot(222); hold on

%ts=y(120,100,4,:); ts=ts(:); plot(ts);
%ts=y(150,106,4,:); ts=ts(:); plot(ts);
%ts=y(250,50,4,:); ts=ts(:); plot(ts);

% agar
%ts=y(125,236,6,:); ts=ts(:); plot(ts);
%ts=y(120,236,8,:); ts=ts(:); plot(ts);
%ts=y(124,112,6,:); ts=ts(:); plot(ts);
%ts=y(123,205,6,:); ts=ts(:); plot(ts);
%%
ye=y(:,:,:,end);
figure
for s=1:size(ye,3)
    hs(s)=subplot(3,4,s);
    img=squeeze(ye(:,:,s));
    imagesc(img);
    %colormap jet
    if s==size(ye,3); hcb=colorbar('east'); end
end
cbpos=get(hcb,'Position');
set(hcb,'Position',[cbpos(1)+0.075 cbpos(2) cbpos(3) cbpos(4)]);

vox=[125,236,6; 120,236,8; 123,205,6];
hs(s+1)=subplot(3,1,3); hold on
for v=1:size(vox,1)
    ts=y(vox(v,1),vox(v,2),vox(v,3),:); ts=ts(:); plot(ts);
% ts=y(125,236,6,:); ts=ts(:); plot(ts);
% ts=y(120,236,8,:); ts=ts(:); plot(ts);
% %ts=y(124,112,6,:); ts=ts(:); plot(ts);
% ts=y(123,205,6,:); ts=ts(:); plot(ts);
end