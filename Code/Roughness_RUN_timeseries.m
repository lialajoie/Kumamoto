% Lia Lajoie
% CSM

% Calculate rouchnesses and make plots
clear all

%% USER INPUTS
PSD_xaxis = 'time'; % "frequency" 

% FREAD_SAC(filename)
% filename = 'HK.HKPS..BHE.M.2016.106.162816.sac';
filename = 'HK.HKPS..BHN.M.2016.106.162816.sac';
% filename = 'HK.HKPS..BHZ.M.2016.106.162816.sac';

filepath = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/Seismic/';
file = [filepath,filename];

[time, data, hdr] = fread_sac(file);
% >> sac_mat = FREAD_SAC(filename)

% limit time
tmin = 62; % sec
tmax = 1000; %3000; % sec
tfind = time >= tmin & time <= tmax;
t = time(tfind);
d = data(tfind);

figure(1)
clf
plot(t,d)

%% TEST PLOT
    figure(1)
    clf
    plot(t,d,'b-')
        title('time series')
        xlabel('time (s)')
        ylabel('offset (units?)')
        set(gca,'fontsize',12,'fontweight','bold')
        grid on

%% CALCULATE PARAMS
prof_space_s = t(2)-t(1);
Si = prof_space_s; %/1000; %/1000; % sampling interval  = profile spacing in km 
xmin = (prof_space_s*2);

%% PSD - Thompson

N = length(d);
x = d;
fs = 1/Si; % sampling frequency
[pxx,f] = pmtm(d,[],length(d),fs); % multitaper PSD
freq = 0:fs/length(x):fs/2;

% linear fit to PSD in log space
    if strcmpi(PSD_xaxis,'time') == 1; % "wavelength" or "frequency"
        % x0_hurst = log10(1./freq);
        x0_hurst = log10(1./f)';
    elseif strcmpi(PSD_xaxis,'frequency') == 1; 
        % x0_hurst = log10(freq);
        x0_hurst = log10(f)';
    end
        y0_hurst = log10(pxx);

    ROLL1 = 1000000000;
    [Hurst_Thom,Roughness_Thom,RSq_Thom,slope_Thom,x_Thom,yCalc_Thom]...
        = fit_lin2log(x0_hurst,y0_hurst,ROLL1);
    Hurst_Thom
    Roughness_Thom
    RSq_Thom
    
% FIGURE 3 - PSD THOMSON
figure (3)
clf
    xt = 10^2;
    yt = 10^0;
    if strcmpi(PSD_xaxis,'time') == 1; % "wavelength" or "frequency"
        loglog(1./freq,(pxx))
        hold on
        xlabel('time (s)')
        set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')
    elseif strcmpi(PSD_xaxis,'frequency') == 1;
        loglog(freq,(pxx))
        hold on
        xlabel('frequency (Hz)')
        set(gca,'fontsize',12,'fontweight','bold')
    end    
        loglog(10.^x_Thom,10.^yCalc_Thom,'r-','linewidth',2)
% 
        text_1 = ['Roughness = ',num2str(Roughness_Thom)];
        text_2 = ['Hurst = ',num2str(Hurst_Thom)];
        text_3 = ['R-squared = ',num2str(RSq_Thom)];
        text(xt,yt,{text_1,text_2,text_3},'fontsize',12,'fontweight','bold')

        ylabel('power spectral density (m^3)')
%         title({'power spectral density - Thomson Multitaper';text_1;text_2}) 
        title({'power spectral density - Thomson Multitaper'})
        set(gca,'fontsize',12,'fontweight','bold')
        legend('PSD','best fit')
        grid on
        
%         ({'line1', 'line2','line3'},)

%% PSD - FFT

% THIBAULT CODE

    %Detrending step:
    profile = d;
    profile = profile - mean(profile); %%
    profile = detrend(profile);
    
    %Tapering step:
    percent = 0.03;
    profile_taper = Ftapering_copy(profile,percent);

% EMILY code

    % z = profile detrended and tapered
    z = profile_taper;
    N = length(z);
    dx = prof_space_s;
    y= fft(z);
    
    % power
    p= y.*conj(y)./(N*dx); 
    
    % put back in dx
    p = p.*dx*dx;
    f = (0:N-1)'/(dx*N);
    p = p(3:N/2);
    f = f(3:N/2);
    
% linear fit to PSD in log space
    if strcmpi(PSD_xaxis,'time') == 1; % "wavelength" or "frequency"
        x0_fft = log10(1./f)';
    elseif strcmpi(PSD_xaxis,'frequency') == 1; 
        x0_fft = log10(f)';
    end
        y0_fft = log10(p);

    ROLL1 = 1000000000;
    ROLL2 = 0.002304;
    [Hurst_fft,Roughness_fft,RSq_fft,slope_fft,x_fft,yCalc_fft]...
        = fit_lin2log(x0_fft,y0_fft,ROLL1);
    Hurst_fft
    Roughness_fft
    RSq_fft
    
% FIGURE 5 - PSD FFT, EMILY
    xt = 10^2;
    yt = 10^0;
figure(5)
    clf
    if strcmpi(PSD_xaxis,'time') == 1; % "wavelength" or "frequency"
        loglog(1./f,p);
        hold on
        xlabel('time (s)')
        set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')
    elseif strcmpi(PSD_xaxis,'frequency') == 1;
        loglog(f,p);
        hold on
        xlabel('frequency (Hz)')
        set(gca,'fontsize',12,'fontweight','bold')
    end
        loglog(10.^x_fft,10.^yCalc_fft,'r-','linewidth',2)
        
        text_1 = ['Roughness = ',num2str(Roughness_fft)];
        text_2 = ['Hurst = ',num2str(Hurst_fft)];
        text_3 = ['R-squared = ',num2str(RSq_fft)];
        text(xt,yt,{text_1,text_2,text_3},'fontsize',12,'fontweight','bold')
        
        ylabel('power spectral density (m^3)')
        title('power spectral density - FFT') 
        set(gca,'fontsize',12,'fontweight','bold')
        legend('PSD','best fit')
        grid on
       

%% WELCH PSD
 % find the step size of X
    [Pyy,Fw] = pwelch(d,[],[],[],fs,'onesided','PSD'); % compute the PSD
  
% linear fit to PSD in log space
   if strcmpi(PSD_xaxis,'time') == 1; % "wavelength" or "frequency"
        x0_welch = log10(1./Fw)';
    elseif strcmpi(PSD_xaxis,'frequency') == 1; 
        x0_welch = log10(Fw)';
    end
        y0_welch = log10(Pyy);

    ROLL1 = 1000000000;
    ROLL2 = 0.002304;
    [Hurst_welch,Roughness_welch,RSq_welch,slope_welch,x_welch,yCalc_welch]...
        = fit_lin2log(x0_welch,y0_welch,ROLL1);
    Hurst_welch
    Roughness_welch
    RSq_welch
    
% FIGURE 4 - PSD WELCH
%     figure(4)
%     plot(log10(F),log10(Pyy)), xlabel('log [\mum^{-1}]'), ylabel('log[PSD[\mum^2/\mum^{-1}]]')
%     title ('Welch PSD')
%     axis tight, grid on
    
    xt = 10^2;
    yt = 10^6;
figure(4)
    clf
    if strcmpi(PSD_xaxis,'time') == 1; % "wavelength" or "frequency"
        loglog(1./Fw,Pyy);
        hold on
        xlabel('time (s)')
        set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')
    elseif strcmpi(PSD_xaxis,'frequency') == 1;
        loglog(Fw,Pyy);
        hold on
        xlabel('frequency (Hz)')
        set(gca,'fontsize',12,'fontweight','bold')        
    end
        loglog(10.^x_welch,10.^yCalc_welch,'r-','linewidth',2)
        
        text_1 = ['Roughness = ',num2str(Roughness_welch)];
        text_2 = ['Hurst = ',num2str(Hurst_welch)];
        text_3 = ['R-squared = ',num2str(RSq_welch)];
        text(xt,yt,{text_1,text_2,text_3},'fontsize',12,'fontweight','bold')
        
        ylabel('power spectral density (m^3)')
        title('power spectral density - Welch') 
        legend('PSD','best fit')
        grid on


% % figure 25
%     Figures25_26_ProfileStatsResults_Kum