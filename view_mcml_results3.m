clear all; close all; clc

load ../data/mcml/Arz

x=[r(1) r(end)];
y=[-z(end) z(end)];

Adraw=cat(2,Abs(end:-1:1,:).',Abs.');

figure;
imagesc(y,x,Adraw);
%surf(Abs);
ylabel('depth (cm)')
xlabel('radius (cm)');
xlim([-2 2]); % 2 cm radial
ylim([0 3]); % 3 cm depth
title('Absorption (J/cm^2)','FontWeight','normal','FontSize',16);
hcb=colorbar;
colormap(jet);
print -dpng ../figures/mcml

%%
% % determine the ROI here
% % look at each depth separately
% thresh=0.99;
% nDepths=size(Abs,2); % 4550
% rc=zeros(nDepths,1);
% for d=1:nDepths
%     thisAbs=Abs(:,d);
%     cs=cumsum(thisAbs);
%     ncs=cs/max(cs);
%     rc(d)=r(find(ncs>thresh,1));
% end
% figure; plot(z,rc);

%%
% find smallest rectangle containing proportion of total absorption
thresh=0.99;
nRads=size(Abs,1); % 2550
nDepths=size(Abs,2); % 4550
rc=zeros(nDepths,1);
totalAbs=sum(sum(Abs));
for d=1:nDepths
        tAbs=Abs;
        tAbs(:,d+1:end)=0; 
        tmp=sum(tAbs,2); % integral over all depths
        tind=find(cumsum(tmp)>thresh*totalAbs,1);
        if isempty(tind)
            minr(d)=NaN;
        else
            minr(d)=r(tind);
        end
end

%%
figure; 
subplot(221); plot(z,minr);
subplot(222); plot(z,minr.*z);
[~,minind]=min(minr.*z)
z(minind)
minr(minind)

% %%
% % here we integrate across the other dimension and arrive at a threshold
% % for z and r separately
% figure;
% subplot(221)
% cs=cumsum(sum(Abs,1));
% ncs=cs/max(cs);
% plot(z, ncs );
% zo=z(find(ncs>thresh,1));
% 
% subplot(222)
% cs=cumsum(sum(Abs,2));
% ncs=cs/max(cs);
% plot(r, ncs );
% ro=r(find(ncs>thresh,1));
% 
% [zo,ro]