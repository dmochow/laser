function [xyzOut,voxOut] = orig2tal(voxIn,Xat,InfoOrig,InfoTlrc)
%Take a point in voxel space of orig BRIK and output xyz talairach

% from voxel orig to XYZ orig
Xout=zeros(3,1);
Xout(1)=InfoOrig.ORIGIN(1)+voxIn(1)*InfoOrig.DELTA(1);
Xout(2)=InfoOrig.ORIGIN(2)+voxIn(2)*InfoOrig.DELTA(2);
Xout(3)=InfoOrig.ORIGIN(3)+voxIn(3)*InfoOrig.DELTA(3);
Xout=[Xout(:);1];

% from XYZ orig to XYZ talairach
xyzOut=inv(Xat)*Xout;
xyzOut=xyzOut(1:3);
xyzOut=xyzOut(:);


% from XYZ talairach to voxel talairach
voxOut=zeros(3,1);
voxOut(1)= (xyzOut(1)-InfoTlrc.ORIGIN(1))/InfoTlrc.DELTA(1);
voxOut(2)= (xyzOut(2)-InfoTlrc.ORIGIN(2))/InfoTlrc.DELTA(2);
voxOut(3)= (xyzOut(3)-InfoTlrc.ORIGIN(3))/InfoTlrc.DELTA(3);
end

