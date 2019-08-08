function deleteRoiFiles(pathToData)
origPath=pwd;
cd(pathToData);
delete roi+tlrc.BRIK;
delete roi+tlrc.HEAD;
delete controlRoi+tlrc.BRIK;
delete controlRoi+tlrc.HEAD;
delete resampled_roi+tlrc.BRIK;
delete resampled_roi+tlrc.HEAD;
delete resampled_control_roi+tlrc.BRIK;
delete resampled_control_roi+tlrc.HEAD;
cd(origPath);
end