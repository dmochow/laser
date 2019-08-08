function arMagSpec = getARspectrum(coeffs,nfft)
if nargin<2, nfft=32; end
coeffs=-coeffs; %
coeffs=cat(1,ones(1,size(coeffs,2),size(coeffs,3)),coeffs);
A=fft(coeffs,nfft,1);
C=1./(A+eps); %invert to get the AR spectrum
arMagSpec=abs(C);