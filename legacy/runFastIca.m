clear all; close all; clc
addpath(genpath('.'))

%% define subject
subjNum=2;

%% load data
BrikFilename=['../data/S0' num2str(subjNum) '/NII/mebold2go/zcat_smooth_deoblique_bold+orig.BRIK'];
Opt.Format = 'matrix';
bold = BrikLoad (BrikFilename, Opt);
[nX,nY,nZ,nT]=size(bold);
%%
MaskFilename=['../data/S0' num2str(subjNum) '/NII/mebold2go/mask+orig.BRIK'];
mask = BrikLoad (MaskFilename, Opt);

%% separate bold into echoes
bold_e1=bold(:,:,1:nZ/3,:);
bold_e2=bold(:,:,nZ/3+1:2*nZ/3,:);
bold_e3=bold(:,:,2*nZ/3+1:nZ,:);


%% mask out non-brain
bold_2D_e1=reshape(bold_e1,[nX*nY*nZ/3 , nT]);
bold_2D_e1=bold_2D_e1(mask(:)>0,:);

bold_2D_e2=reshape(bold_e2,[nX*nY*nZ/3 , nT]);
bold_2D_e2=bold_2D_e2(mask(:)>0,:);

bold_2D_e3=reshape(bold_e3,[nX*nY*nZ/3 , nT]);
bold_2D_e3=bold_2D_e3(mask(:)>0,:);

%%
bold_2D=cat(1,bold_2D_e1,bold_2D_e2);
bold_2D=cat(1,bold_2D,bold_2D_e3);

%% demean and variance normalize
bold_2D_z = zscore(bold_2D,0,2); % normalize along the time dimension
clear bold;
clear bold_e1;
clear bold_e2;
clear bold_e3;


% %% estimate dimensionality of data
% % standard PCA
% R=bold_2D_z'*bold_2D_z;
% [V,D]=eig(R);
% %%
% nBrainVoxels=size(bold_2D_z,1)/3;
% nComp=size(V,2);
% 
% for c=1:nComp
%     [c,v]
%     W=V(:,c);
%     comp = bold_2D_z*W;
%     
%     comp_e1=comp(1:nBrainVoxels,:);
%     comp_e2=comp(nBrainVoxels+1:2*nBrainVoxels,:);
%     comp_e3=comp(2*nBrainVoxels+1:3*nBrainVoxels,:);
%    
%     for v=1:nBrainVoxels
%         y=[comp_e1(v);comp_e2(v);comp_e3(v)];
%         X=[ (1:3)', ones(3,1)];
%         [b,bint,r,rint,stats]=regress(y,X);
%         error=y-X*b;
%         error_0=y;
%         fval(v)=sum(error.^2)/sum(error_0.^2);
%         %fvals(v)=stats(2);
%         %X=[ [0;0;0], ones(3,1)];
%         %[b,bint,r,rint,stats]=regress(y,X);
%         %kappa(v)=stats(2);
%     end
%     
%     kappa(c)=(fval* (comp_e1.^2))/sum(comp_e1.^2);
%     
% end


%% estimate dimensionality of data
%K=size(bold_z_2D,2)-1;
% K=20;
% [coeff,score,pcvar] = ppca(bold_2D_z,K);


% 
% %%
% %figure;
% %semilogy(diag(D));
% 

% 
% %% run ICA -- this shit is crashing hard
% %[bold_ica, A, W] = fastica (bold_z_2D.');
