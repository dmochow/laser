clear all; close all; clc
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));

% load first frame from functional
% NB: suBOLD_AP.nii is already registered to the anatomy via the header
boldFilename='/Users/jacek/Documents/MATLAB/LLLT/data/S01/NII/uBOLD_AP.nii';
boldNii = load_untouch_nii(boldFilename,1);
boldNiiFull = load_untouch_nii(boldFilename);
[bold_sx, bold_sy, bold_sz]=size(boldNii.img);
bold_ts=double(boldNiiFull.img);
vox2ras_bold=[boldNii.hdr.hist.srow_x; boldNii.hdr.hist.srow_y; boldNii.hdr.hist.srow_z];
vox2ras_bold=cat(1,vox2ras_bold,[0 0 0 1]);

% data from unwrap and realign
load('/Users/jacek/Documents/MATLAB/LLLT/data/S01/NII/BOLD_AP_uw.mat','ds');

%%

% load anatomical (i.e., grey matter mask)
greyMaskFilename='/Users/jacek/Documents/MATLAB/LLLT/data/S01/NII/c3T1w_MPR_BIC_v1.nii';
greyNii = load_untouch_nii(greyMaskFilename);
greyMask=double(greyNii.img);
vox2ras_T1=[greyNii.hdr.hist.srow_x; greyNii.hdr.hist.srow_y; greyNii.hdr.hist.srow_z];
vox2ras_T1=cat(1,vox2ras_T1,[0 0 0 1]);

% get sizes of volumens
bold2anat=zeros(bold_sx,bold_sy,bold_sz);
%[anat_sx, anat_sy, anat_sz]=size(greyNii.img);
anat_sz=size(greyNii.img);

%% map each voxel in bold to corresponding voxel in anatomy
for ix=1:bold_sx
    for iy=1:bold_sy
        for iz=1:bold_sz
            
            thisVoxel=[ix;iy;iz;1];
            bold_ras=vox2ras_bold*thisVoxel;
            anat_vox=inv(vox2ras_T1)*bold_ras;
            anat_vox=round(anat_vox);
            try
                thisInd = sub2ind(anat_sz,anat_vox(1),anat_vox(2),anat_vox(3));
            catch
                thisInd=NaN;
            end
            bold2anat(ix,iy,iz)=thisInd;
            
        end
    end
end


%% figure out for each bold voxel whether it's grey matter
for ix=1:bold_sx
    for iy=1:bold_sy
        for iz=1:bold_sz
            
            thisGreyInd=bold2anat(ix,iy,iz);
            if ~isnan(thisGreyInd)
                boldIsGrey(ix,iy,iz)=greyNii.img(thisGreyInd)>0;
            else
                boldIsGrey(ix,iy,iz)=0;
            end
            
        end
    end
end

%% now grab time series from all grey voxels
tmp=boldIsGrey(:);
bold_ts_2D=reshape(bold_ts,[bold_sx*bold_sy*bold_sz size(bold_ts,4)]);
grey_ts_2D=bold_ts_2D(tmp==1,:);

mu_grey_ts_2D=mean(grey_ts_2D,1);

%%
% regress out deformation fields
grey_ts_2D_out = regressOut(grey_ts_2D,ds.q.');
mu_grey_ts_2D_out=mean(grey_ts_2D_out);
%%
figure;
subplot(211);
plot(mu_grey_ts_2D);
subplot(212);
plot(mu_grey_ts_2D_out);

% %% run anova on grey voxels
% nGreyVoxels=size(grey_ts_2D,1);
% g1=1:1200;
% g2=cell(1,1200);
% for i=1:1200
%     if i<401
%         g2{i}='off';
%     elseif i<801
%         g2{i}='on';
%     else
%         g2{i}='post';
%     end
% end
% 
% %%
% for vox=2000%nGreyVoxels
%     y=grey_ts_2D(vox,:)';
%     [p,tbl]=anovan(y,{g1,g2});
% end
    