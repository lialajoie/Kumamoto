% Lia Lajoie
% CSM
% 12 April 2018

clear all

% plot Kumamoto waveform data and calculate roughness.

% FREAD_SAC(filename)
% filename = 'HK.HKPS..BHE.M.2016.106.162816.sac';
filename = 'HK.HKPS..BHN.M.2016.106.162816.sac';
% filename = 'HK.HKPS..BHZ.M.2016.106.162816.sac';

filepath = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/Seismic/';
file = [filepath,filename];

[t, data, hdr] = fread_sac(file);
% >> sac_mat = FREAD_SAC(filename)

figure(1)
clf
plot(t,data)