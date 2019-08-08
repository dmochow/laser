% simulate the bold signal
clear all; close all; clc

t=0:0.001:0.1; % time
t2s=0.05; % T2*
So=1;
delSo=0.05;
delT2s=0.03;

s=So*exp(-t/t2s);
s1=(So+delSo)*exp(-t/t2s);
s2=So*exp(-t/(t2s*(1+delT2s)));
s3=(So+delSo)*exp(-t/(t2s*(1-delT2s)));

dels1=s1-s;
dels2=s2-s;
dels3=s3-s;

figure;
% subplot(231); hold on
% plot(t,s,'k');
% plot(t,s1,'r');
% plot(t,s2,'g');
subplot(231); hold on
plot(t,dels1,'k');
subplot(232); hold on
plot(t,dels2,'k');
subplot(233); hold on
plot(t,dels3,'k');

subplot(234); hold on
plot(t,dels1./s1,'k');
subplot(235); hold on
plot(t,dels2./s2,'k');
subplot(236); hold on
plot(t,dels3./s3,'k');