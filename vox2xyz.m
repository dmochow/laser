function voxOut = xyz2vox(xyzIn,Info)
%convert from xyz to voxel
voxOut=zeros(3,1);
voxOut(1)= (xyzIn(1)-Info.ORIGIN(1))/Info.DELTA(1);
voxOut(2)= (xyzIn(2)-Info.ORIGIN(2))/Info.DELTA(2);
voxOut(3)= (xyzIn(3)-Info.ORIGIN(3))/Info.DELTA(3);
end

