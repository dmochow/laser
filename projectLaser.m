function roiMask = projectLaser(brainMaskFilename,laserOrigin,urf,ulf,urb,lrf,radMask,depthMask)

if nargin<7
radMask=19.2; % lateral laser reach in mm
depthMask=26.5; % depth laser reach in mm
end
% project the real laser into the mri
[err, mask, Info, ErrMessage] = BrikLoad (brainMaskFilename);
mridim=size(mask);
[x,y,z]=ndgrid(0:mridim(1)-1,0:mridim(2)-1,0:mridim(3)-1);
LOx=laserOrigin(1)*ones(size(x));
LOy=laserOrigin(2)*ones(size(y));
LOz=laserOrigin(3)*ones(size(z));

%%
% translation
T=eye(4);
T(1,4)=-laserOrigin(1);
T(2,4)=-laserOrigin(2);
T(3,4)=-laserOrigin(3);

xt=T(1,1)*x+T(1,2)*y+T(1,3)*z+T(1,4)*1;
yt=T(2,1)*x+T(2,2)*y+T(2,3)*z+T(2,4)*1;
zt=T(3,1)*x+T(3,2)*y+T(3,3)*z+T(3,4)*1;

%%
% rotation
u=(urf-ulf)/norm(urf-ulf);  % LR
v=(urf-urb)/norm(urf-urb);  % PA
w=(urf-lrf)/norm(urf-lrf);  % IS
R=eye(4);
R(1:3,1)=u;
R(1:3,2)=v;
R(1:3,3)=w;
R=inv(R);

xr=R(1,1)*xt+R(1,2)*yt+R(1,3)*zt+R(1,4)*1;
yr=R(2,1)*xt+R(2,2)*yt+R(2,3)*zt+R(2,4)*1;
zr=R(3,1)*xt+R(3,2)*yt+R(3,3)*zt+R(3,4)*1;

%%
roiMask=zeros(size(mask));
roiMask( (xr.^2+zr.^2)<radMask.^2 & yr>0 & yr<depthMask )=1;

%
% r = sqrt(xr.^2+zr.^2)
% z= yr