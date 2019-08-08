clear all; close all; clc

[X,Y,Z]=meshgrid(1:100,1:100,1:100);

% define 8 points
llf=[80 80 80];
ulf=[80 80 86];
llb=[80 86 80];
ulb=[80 86 86];

lrf=[86 80 80];
urf=[86 80 86];
lrb=[86 86 80];
urb=[86 86 86];

alpha=15*pi/180;
% R=[    1.0000         0         0
%          0    cos(alpha)   -sin(alpha)
%          0    sin(alpha)    cos(alpha) ];
% 
% llb=(R*llb.').';
% llf=(R*llf.').';
% ulf=(R*ulf.').';
% ulb=(R*ulb.').';
% lrf=(R*lrf.').';
% urf=(R*urf.').';
% lrb=(R*lrb.').';
% urb=(R*urb.').';

xslice = 80:90; 
yslice = 80:90; 
zslice = 80:90;

% define "left" plane
cross_l=cross((llf-llb),(llf-ulf));
xo=llf;
X0=repmat(xo(1),[size(X,1) size(X,2) size(X,3)]);
Y0=repmat(xo(2),[size(Y,1) size(Y,2) size(Y,3)]);
Z0=repmat(xo(3),[size(Z,1) size(Z,2) size(Z,3)]);
isRight= ( cross_l(1)*(X-X0)+cross_l(2)*(Y-Y0)+cross_l(3)*(Z-Z0) ) > 0;

% define "right" plane
cross_r=cross((lrf-lrb),(lrf-urf));
xo=lrf;
X0=repmat(xo(1),[size(X,1) size(X,2) size(X,3)]);
Y0=repmat(xo(2),[size(Y,1) size(Y,2) size(Y,3)]);
Z0=repmat(xo(3),[size(Z,1) size(Z,2) size(Z,3)]);
isLeft= ( cross_r(1)*(X-X0)+cross_r(2)*(Y-Y0)+cross_r(3)*(Z-Z0) ) < 0;

% define "front" plane
cross_f=cross((llf-ulf),(llf-lrf));
xo=llf;
X0=repmat(xo(1),[size(X,1) size(X,2) size(X,3)]);
Y0=repmat(xo(2),[size(Y,1) size(Y,2) size(Y,3)]);
Z0=repmat(xo(3),[size(Z,1) size(Z,2) size(Z,3)]);
isFront= ( cross_f(1)*(X-X0)+cross_f(2)*(Y-Y0)+cross_f(3)*(Z-Z0) ) < 0;


% define "upper" plane
cross_u=cross((ulf-ulb),(ulf-urf));
xo=ulf;
X0=repmat(xo(1),[size(X,1) size(X,2) size(X,3)]);
Y0=repmat(xo(2),[size(Y,1) size(Y,2) size(Y,3)]);
Z0=repmat(xo(3),[size(Z,1) size(Z,2) size(Z,3)]);
isBelow= ( cross_u(1)*(X-X0)+cross_u(2)*(Y-Y0)+cross_u(3)*(Z-Z0) ) > 0;

% define "lower" plane
cross_l=cross((llf-llb),(llf-lrf));
xo=llf;
X0=repmat(xo(1),[size(X,1) size(X,2) size(X,3)]);
Y0=repmat(xo(2),[size(Y,1) size(Y,2) size(Y,3)]);
Z0=repmat(xo(3),[size(Z,1) size(Z,2) size(Z,3)]);
isAbove= ( cross_l(1)*(X-X0)+cross_l(2)*(Y-Y0)+cross_l(3)*(Z-Z0) ) < 0;



isAll=and(and(and(and(isLeft,isRight),isFront),isBelow),isAbove);


% figure
% subplot(221);
% slice(X,Y,Z,double(isRight),xslice,yslice,zslice)
% colormap bone
% 
% subplot(222);
% slice(X,Y,Z,double(isLeft),xslice,yslice,zslice)
% colormap bone
% 
% subplot(223);
% slice(X,Y,Z,double(isFront),xslice,yslice,zslice)
% colormap bone
% 
% subplot(224);
% slice(X,Y,Z,double(isAll),xslice,yslice,zslice)
% colormap bone


%%
figure; hold on
plot3(llf(1),llf(2),llf(3),'*k');
plot3(ulf(1),ulf(2),ulf(3),'*k');
plot3(llb(1),llb(2),llb(3),'*k');
plot3(ulb(1),ulb(2),ulb(3),'*k');
plot3(lrf(1),lrf(2),lrf(3),'*k');
plot3(urf(1),urf(2),urf(3),'*k');
plot3(lrb(1),lrb(2),lrb(3),'*k');
plot3(urb(1),urb(2),urb(3),'*k');
xlim([1 100]); ylim([1 100]); zlim([1 100]);
%slice(X,Y,Z,double(isAll),xslice,yslice,zslice)
H = vol3d('CData',double(isAll));
colormap(jet);

% figure; 
% hold on


