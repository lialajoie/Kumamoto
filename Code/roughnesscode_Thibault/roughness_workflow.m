clear all
%close all

% load toto.mat

% SLIP FROM THIBAULT ET AL
    % slip_EQ = 'Borah_Peak_XZ.mat';
    % slip_EQ = 'Hector_mine_XZ.mat';
    slip_EQ = 'Landers_XZ.mat';
    % slip_EQ = 'Kumamoto_XZ.mat';
    % slip_EQ = 'Kumamoto_LLajoie_40m_XZ.mat';
    % slip_EQ = 'Landers_disp_data_no_header.csv';
    % slip_EQ = 'Superstition_Hills_XZ.mat';
    slip_dir = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/';
    slip_file = [slip_dir,slip_EQ];

    if strcmpi(slip_EQ,'Landers_XZ.mat') == 1
        S = open(slip_file);
        Z = S.Z';
        X = S.X';
    elseif strcmpi(slip_EQ,'Kumamoto_XZ.mat') == 1
        S = open(slip_file);
        X = S.x';
        Z = S.z';
    elseif strcmpi(slip_EQ,'Kumamoto_LLajoie_40m_XZ.mat') == 1
        S = open(slip_file);
        X = S.distalong_transect_km*1000;
        Z = S.disp_m;  
        
%         % pseudo- full-slip
%         add_dist = distalong_transect_m + max(distalong_transect_m) + (distalong_transect_m(2)-distalong_transect_m(1));
%         distalong_transect_m = [distalong_transect_m add_dist];
%         disp_m = [disp_m fliplr(disp_m)]
 
    elseif strcmpi(slip_EQ,'Landers_disp_data_no_header.csv') == 1
        S = importdata(slip_file);
        distalong_transect_km = S(:,1);
        X = distalong_transect_km*1000';
        Z = S(:,2)';
    else
        S = open(slip_file);
        Z = S.Z;
        X = S.X;
    end
    
% L_data = load('Landers_disp_data_no_header');
% distalong_transect_km = L_data(:,1);
% X = distalong_transect_km*1000';
% Z = L_data(:,2)';

figure(1)
clf
plot(X,Z)

dx = X(2) - X(1);
[nz,nx]=size(Z);

for i=1:nz;
profile=Z(i,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Detrending step:
profile = profile - mean(profile); %%
profile = detrend(profile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tapering step:
percent = 0.03;
profile_taper = Ftapering(profile,percent);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fourier
if i==1
[k3,ft3] = F_ps_profile1D_new(profile,dx);
ave_ps   = ft3;

else
[k3,ft3] = F_ps_profile1D_new(profile,dx);
ave_ps   = ave_ps + ft3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

clear ft3
ft3 = ave_ps/nz;

k3_final = k3;
ft3_final = ft3;

%I apply the re-sampling
%minK and maxK can be varied if we want to keep only the frequencies lower than those affected by the acquisition noise
% minK = k3_final(2);
% maxK = k3_final(end);

%% LIA ADD
%numbers from (Candela et al., 2011)
minK = 10^-5;
maxK = 4e-04;

% plot lims 
maxlim = log10(1/10^-5)
minlim = log10(1/4e-04)
%% 

nbins = 20;
[tX, yy1] = Ftrendfft2(k3_final,ft3_final,minK,maxK,nbins);

Wavelength = log10(1./tX)';
PSD = log10(yy1);

% Wavelength = 1./tX;
% PSD = yy1;

figure(2)
clf
% loglog(Wavelength,PSD,'co-');hold on
plot(Wavelength,PSD,'co-');hold on

%%
