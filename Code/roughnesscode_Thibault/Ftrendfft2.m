function [tX, yy1] = Ftrendfft2(k1,fft1,minK,maxK,nbins)
%function [tX, tY] = Ftrendfft(k1,fft1,minK,maxK,nbins)
% Calcule un trend en loi puissance sur une fft, avec binning
nbins = nbins +1;
bins    = log10(minK):(log10(maxK) -log10(minK))/(nbins-1):log10(maxK);
k_bin   = zeros(1,(nbins-1));
fft_bin = zeros(1,(nbins-1));
for i=1:(nbins-1)
    J          = find ((log10(k1) >= bins(i)) & (log10(k1) < bins(i+1)));
    k_bin(i)   = mean(k1(J));
    fft_bin(i) = mean(fft1(J));
end
isnan_fft = isnan(fft_bin);
J         = find(isnan_fft == 1);
k_bin(J)  = [];
fft_bin(J)= [];

tX        = log10(k_bin);
yy1       = log10(fft_bin);
p1        = polyfit(tX,yy1,1);
tY        = 10^p1(2)*k_bin.^p1(1);
tX        = 10.^tX;
yy1       = 10.^yy1;
H = -(1 +p1(1))/2
D = (5 + p1(1))/2 % Lia addition
