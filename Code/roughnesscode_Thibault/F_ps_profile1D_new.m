function [k,ps] = F_ps_profile1D(profX,pixel_size)
% compute the power spectrum of the profile

nx   = length(profX);
fftx = fft(profX);
ps   = fftx .* conj(fftx)/length(profX);
k    = (1/pixel_size)./2*linspace(0,1,(floor((length(profX))/2)+1));
ps   = ps(1:(length(profX)/2+1));
ps   = ps.*(pixel_size);







