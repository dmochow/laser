% 07/20/18
% 
% find smallest rectangle containing proportion of total absorption
%
clear all; close all; clc
thresh=0.95;
load ../data/mcml/Arz

nRads=size(Abs,1); % 2550
nDepths=size(Abs,2); % 4550
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
% figure; 
% subplot(221); plot(z,minr);
% subplot(222); plot(z,minr.*z);
[~,minind]=min(pi*minr.*minr.*z)
zo=z(minind);
ro=minr(minind);

%%
fid=fopen(['../data/minAreaROI_' num2str(thresh*100) 'percent'],'w');
fwrite(fid,['depth = ' num2str(zo) newline]);
fwrite(fid,['radius = ' num2str(ro) newline]);
fclose(fid);

