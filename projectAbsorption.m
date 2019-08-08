function Amri = projectAbsorption(brainMaskFilename,laserOrigin,urf,ulf,urb,lrf)

%%%% LOAD ABSORPTION
% load('../data/precomputed/RZA.mat','R','Z','A');
% Ra=R; Za=Z;
% ra=Ra(1,:).'; za=Za(:,1);

load('../data/precomputed/RZAfull_092418.mat','R','Z','A');
ra=R(:); za=Z(:);

if nargin<6, error('need 6 arguments'); end
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
% now index from xr,yr,zr to ra,za
% define query points
rq=sqrt(xr(:).^2+zr(:).^2); 
zq=yr(:);
Aq=zeros(size(zq)); % absorption on mri
nQ=size(zq,1);
mapR=NaN(nQ,1);
mapZ=NaN(nQ,1);

% define limits
maxZ=max(za);
maxR=max(ra);

for q=1:nQ
    if zq(q)>0 && zq(q)<maxZ && rq(q)<maxR
        [~,rind]=min(abs(rq(q)-ra));
        mapR(q)=rind;

        [~,zind]=min(abs(zq(q)-za));
        mapZ(q)=zind;

        Aq(q)=A(zind,rind);
    end
end

% reshape Aq
Amri=reshape(Aq,[mridim(1) mridim(2) mridim(3)]);
