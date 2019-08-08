clear all; close all; clc

a=[1 -0.2 0.4];
[H,F] = freqz(1,a,4);

A=fft(a,8);
A=A(:)
B=1./A;

figure;
subplot(211);
plot(abs(H));
subplot(212);
plot(abs(B));
xlim([1 4])
[H B(1:4)]