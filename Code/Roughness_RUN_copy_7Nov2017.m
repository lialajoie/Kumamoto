% Lia Lajoie
% CSM

% Calculate rouchnesses and make plots
clear all

%% USER INPUTS
L_data = load('Landers_disp_data_no_header');
distalong_transect_km = L_data(:,1);
distalong_transect_m = distalong_transect_km*1000;
disp_m = L_data(:,2);

%% TEST PLOT
    figure(1)
    clf
    plot(distalong_transect_km,disp_m)
    title('Slip Dist')

%% CALCULATE PARAMS
prof_space_m = distalong_transect_m(2)-distalong_transect_m(1);
Si = prof_space_m; %/1000; %/1000; % sampling interval  = profile spacing in km 
xmin = (prof_space_m*2);

% %% PSD - Thompson
% cum_Disp = disp_m;
% N = length(cum_Disp);
% x = cum_Disp;
% fs = 1/Si; % sampling frequency
% % xK = 0:1/fs:1-1/fs*(length(cum_Disp)); % vector of x-values
% %[pxx,f,pxxc] = pmtm(cum_Disp,[],length(cum_Disp),fs,'ConfidenceLevel',0.95); % multitaper PSD
% [pxx,f] = pmtm(cum_Disp,[],length(cum_Disp),fs); % multitaper PSD
% [pxxP,fP] = periodogram(cum_Disp,[],length(cum_Disp),fs);
% freq = 0:fs/length(x):fs/2;
% 
% % linear fit to PSD in log space
%     x0 = log10(1./freq);
%     y0 = log10(pxx);
%     
%     ROLL1 = 1000000000;
%     ROLL2 = 1000000000;
%     
%     % limit to wavelengths < fault length
%     xfind_1 = x0 <= log10(ROLL1);
%     x_1 = x0(xfind_1);
%     y_1 = y0(xfind_1);
%     X_1 = [ones(length(x_1),1) x_1'];
% 
%     b_1 = X_1\y_1;
%     b_1(2)
%     yCalc_Thom1 = X_1*b_1;
%     Rsq2_Thom1 = 1 - sum((y_1 - yCalc_Thom1).^2)/sum((y_1 - mean(y_1)).^2)
%     Roughness_Thom1 = (5-b_1(2))/2
%     Hurst_Thom1 = ((b_1(2))-1)/2
%     
%     % limit to wavelengths < rollover (linear section)
%     xlim_2 = [xmin ROLL2];
%     xfind_2 = x0 <= log10(xlim_2(2)) & x0 >= log10(xlim_2(1));
%     x_2 = x0(xfind_2);
%     y_2 = y0(xfind_2);
%     X_2 = [ones(length(x_2),1) x_2'];
% 
%     b_2 = X_2\y_2;
%     b_2(2)
%     yCalc_Thom2 = X_2*b_2;
%     Rsq2_Thom2 = 1 - sum((y_2 - yCalc_Thom2).^2)/sum((y_2 - mean(y_2)).^2)
%     Roughness_Thom2 = (5-b_2(2))/2
%     Hurst_Thom2 = (b_2(2)-1)/2
%     
%     figure (3)
%     clf
%     loglog(1./f,(pxx))
%     hold on
%     loglog(10.^x_1,10.^yCalc_Thom1,'k-','linewidth',2)
% %     hold on
% %     loglog([xlim_1 xlim_1],[10^-2 10^2])
%     hold on
%     loglog(10.^x_2,10.^yCalc_Thom2,'r--','linewidth',2)
%    % xlim([0 30])
%     grid on


%% PSD - FFT
x = cum_Disp;
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
% wave = xK;
wave = 1./freq;
freq = 0:fs/length(x):fs/2;

% linear fit to PSD in log space
    x0 = log10(1./freq);
    y0 = log10(psdx);
    
    ROLL1 = 1000000000; % if you want to fit linear trend to specific parts of curve
    ROLL2 = 1000000000;
    
    % limit to wavelengths < fault length
    xfind_1 = x0 <= log10(ROLL1);
    x_1 = x0(xfind_1);
    y_1 = y0(xfind_1);
    X_1 = [ones(length(x_1),1) x_1'];

    b_1 = X_1\y_1;
    b_1(2)
    yCalc_1 = X_1*b_1;
    Rsq2_1 = 1 - sum((y_1 - yCalc_1).^2)/sum((y_1 - mean(y_1)).^2)
    Roughness_1 = (5-b_1(2))/2
    Hurst_1 = ((b_1(2))-1)/2
    
    % limit to wavelengths < rollover (linear section)
    xlim_2 = [xmin ROLL2];
    xfind_2 = x0 <= log10(xlim_2(2)) & x0 >= log10(xlim_2(1));
    x_2 = x0(xfind_2);
    y_2 = y0(xfind_2);
    X_2 = [ones(length(x_2),1) x_2'];

    b_2 = X_2\y_2;
    b_2(2)
    yCalc_2 = X_2*b_2;
    Rsq2_2 = 1 - sum((y_2 - yCalc_2).^2)/sum((y_2 - mean(y_2)).^2)
    Roughness_2 = (5-b_2(2))/2
    Hurst_2 = (b_2(2)-1)/2

figure(2)
    loglog(wave,(psdx))
    grid on
    title('Periodogram Using FFT')
    xlabel('Wavelength (m)')
    ylabel('Power/Frequency (m^3)')


%% WELCH PSD
 % find the step size of X
    [Pyy,F] = pwelch(cum_Disp,[],[],[],fs,'onesided','PSD'); % compute the PSD
    
    figure(4)
    plot(log10(F),log10(Pyy)), xlabel('log [\mum^{-1}]'), ylabel('log[PSD[\mum^2/\mum^{-1}]]')
    title ('Welch PSD')
    axis tight, grid on
    
%% PLOT

% figure 5 - PSD - FFT
    figure(5)
    clf
    loglog(1./freq,(psdx),'linewidth',2)

    hold on
    loglog(10.^x_1,10.^yCalc_1,'k-','linewidth',2)
%     hold on
%     loglog([xlim_1 xlim_1],[10^-2 10^2])
    hold on
    loglog(10.^x_2,10.^yCalc_2,'k--','linewidth',2)
%     hold on
%     loglog([xlim_2 xlim_2],[10^-2 10^2])
%     title('Multitaper power spectral density - Thompson') 
    xlabel('wavelength(km)')
    ylabel('power spectral density')
    set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')
    xlim([00 70000])
    legtext_1 = ['Roughness = ',num2str(Roughness_1)];
    legtext_2 = ['Roughness = ',num2str(Roughness_2)];
    legend('PSD','best fit',legtext_1,legtext_2)
    pbaspect([2.5 1 1]) 
    grid on
    print('roughness_40','-dpdf')
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print('roughness','-dpng','-r600')
    
%% FIGURE 5
figure(5)
plot(1./f,10*log10(pxx))



figure(7)
plot(log10(f),log10(pxx))
xlim([-1.2 1.2])
grid on
% 
% % figure 25
%     Figures25_26_ProfileStatsResults_Kum