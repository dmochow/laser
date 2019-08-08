function h = colorBrain(anatShow,isSig1,isSig2,brainMask,crop,nColors)
if nargin<6, nColors=64; end
if nargin<5, crop=0; end
anatShow=anatShow.';
isSig1=isSig1.';
isSig2=isSig2.';
brainMask=brainMask.';

anatShow(~brainMask)=NaN;
 
h=pcolor(anatShow); set(h,'EdgeColor','none');
cm=colormap(cat(1,gray(nColors),[255 218 185]/255));
%cm=colormap(cat(1,gray(nColors),[0 1 0]));
cm=colormap(cat(1,cm,[1 0 0]));
cmin = min(anatShow(:));
cmax = max(anatShow(:));
C1 = min(nColors,round((nColors-1)*(anatShow-cmin)/(cmax-cmin))+1); 
newCData=C1;
newCData(isSig1)=nColors+1;
newCData(isSig2)=nColors+2;
set(h,'CData',newCData);
caxis([1 nColors+2]);
axis ij; axis equal; axis off;


if crop
    tmp=mean(~isnan(anatShow));
    leftMargin=find(tmp~=0,1);
    rightMargin=find(tmp~=0,1,'last');
    tmp=mean(~isnan(anatShow),2);
    bottomMargin=find(tmp~=0,1,'last');
    xlim([leftMargin rightMargin]); 
    yl=ylim;
    ylim([yl(1) bottomMargin]);
end
    