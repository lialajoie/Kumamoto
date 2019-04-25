% Lia Lajoie
% CSM

% Calculate rouchnesses and make plots
clear all

%% USER INPUTS
PSD_xaxis = 'wavelength'; % 'frequency' or 'wavelength'
slide_size = .5; % in m - how much windows slide
 
%% SELECT EARTHQUAKE
% LANDERS DATA
% L_data = load('Landers_disp_data_no_header');
% distalong_transect_km = L_data(:,1);
% distalong_transect_m = distalong_transect_km*1000;
% disp_m = L_data(:,2);

% FAULT MORPHOLOGY FROM THIBAULT ET AL
%     EQ = 'Borah_Peak_XZ.mat';
%     EQ = 'Hector_mine_XZ.mat';
%     EQ = 'Kumamoto_XZ.mat';
%     EQ = 'Landers_XZ.mat';
%     EQ = 'Superstition_Hills_XZ.mat';
    
% FAULT SLIP DATASETS
    EQ = 'Kumamoto_LLajoie_40m_XZ.mat';
   % EQ = 'Landers_disp_data_no_header.csv'; % from Milliner
   % EQ = 'HectorMine_Milliner.csv'; % from Milliner
   % EQ = 'Balochistan_Vallage_FarField.txt'; % from Vallage
 
%% RETRIEVE EARTHQUAKE DATA
    EQ_dir = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/';
    EQ_file = [EQ_dir,EQ];
    
    if strcmpi(EQ,'Kumamoto_XZ.mat') == 1
        S = open(EQ_file);
        distalong_transect_m = S.x;
        disp_m = S.z;
        EQ_name = 'Kumamoto - Candela';
    elseif strcmpi(EQ,'Kumamoto_LLajoie_40m_XZ.mat') == 1
        S = open(EQ_file);
        distalong_transect_m = S.distalong_transect_km*1000;
        disp_m = S.disp_m;  
        EQ_name = 'Kumamoto - Lajoie';
        % pseudo- full-slip
%         add_dist = distalong_transect_m + max(distalong_transect_m) + (distalong_transect_m(2)-distalong_transect_m(1));
%         distalong_transect_m = [distalong_transect_m add_dist];
%         disp_m = [disp_m fliplr(disp_m)]
    elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
        S = importdata(EQ_file);
        distalong_transect_km = S(:,1);
        distalong_transect_m = distalong_transect_km*1000;
        disp_m = S(:,2);
        EQ_name = 'Landers - Milliner';
    elseif strcmpi(EQ,'HectorMine_Milliner.csv') == 1
        S = importdata(EQ_file);
        disp_m = S(:,1);
        distalong_transect_m = [0:1:length(disp_m)-1]*100;
        distalong_transect_km = distalong_transect_m/1000;
        EQ_name = 'Hector Mine - Milliner';
    elseif strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        S = importdata(EQ_file);
        distalong_transect_km = S(:,1); % distance from south
        distalong_transect_km = sort(distalong_transect_km);
        distalong_transect_m = distalong_transect_km*1000;
        disp_m = S(:,2); % fault parallel
        EQ_name = 'Balochistan - Vallage'        
    else
        S = open(EQ_file);
        distalong_transect_m = S.X;
        disp_m = S.Z;
        EQ_name = 'Tibo EQ';
    end
    distalong_transect_km = distalong_transect_m/1000;

%% CALCULATE PARAMS
prof_space_m = mean(distalong_transect_m(2:end)-distalong_transect_m(1:end-1));
std_profspace = std(distalong_transect_m(2:end)-distalong_transect_m(1:end-1));
Si = prof_space_m; %/1000; %/1000; % sampling interval  = profile spacing in km 
xmin = (prof_space_m*2);

%% PSD - Thompson

cum_Disp = disp_m;
N = length(cum_Disp);
x = cum_Disp;
fs = 1/Si; % sampling frequency
[pxx,f] = pmtm(cum_Disp,[],length(cum_Disp),fs); % multitaper PSD
freq = 0:fs/length(x):fs/2;

% linear fit to PSD in log space
    if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
        % x0_hurst = log10(1./freq);
        x0_hurst = log10(1./f)';
    elseif strcmpi(PSD_xaxis,'frequency') == 1; 
        % x0_hurst = log10(freq);
        x0_hurst = log10(f)';
    end
        y0_hurst = log10(pxx);

    ROLL1 = 2e10;
%     if strcmpi(fit_type,'whole') == 1
        [H_Thom,D_Thom,RSq_Thom,slope_Thom,x_Thom,yCalc_Thom]...
            = fit_lin2log(x0_hurst,y0_hurst,ROLL1);
%     elseif strcmpi(fit_type,'sliding') == 1
%         x0_hurst_sort = sort(x0_hurst);
%         rm_infinite = isfinite(x0_hurst_sort);
%         x0_hurst_sort = x0_hurst_sort(rm_infinite);
%         win_start = [x0_hurst_sort(1):slide_size*slide_overlap:x0_hurst_sort(end)];
%         win_end = [win_start(2:end) x0_hurst_sort(end)];
%         num_win = length(win_start);
%         for kk = 1:num_win
%             find_win = (x0_hurst >= win_start(kk) & x0_hurst <= win_end(kk));
%             x0_hurst_win = 
%             y0_hurst_win
%             [Hurst_Thom_slide,Roughness_Thom_slide,RSq_Thom_slide,slope_Thom_slide,...
%                 x_Thom_slide,yCalc_Thom_slide] = fit_lin2log(x0_hurst_win,y0_hurst_win,ROLL1);
%         
%         
%         Hurst_Thom_vect = 
%         Roughness_Thom
%         RSq_Thom;
%     end
    
    H_Thom
    D_Thom
    RSq_Thom;
    
% FIGURE 3 - PSD THOMSON
figure (3)
clf
    xt = 10^4.2;
    yt = 10^0;
    if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
        loglog(1./freq,(pxx))
        hold on
        xlabel('wavelength (m)')
        set(gca,'fontsize',12,'fontweight','bold')
    elseif strcmpi(PSD_xaxis,'frequency') == 1;
        loglog(freq,(pxx))
        hold on
        xlabel('wavenumber (m^{-1})')
        set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')
    end    
        loglog(10.^x_Thom,10.^yCalc_Thom,'r-','linewidth',2)
% 
        text_1 = ['Roughness = ',num2str(D_Thom)];
        text_2 = ['Hurst = ',num2str(H_Thom)];
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
    profile = cum_Disp;
    profile = profile - mean(profile); %%
    profile = detrend(profile);
        
    %Tapering step:
    percent = 0.03;
    profile_taper = Ftapering_copy(profile,percent);

% EMILY code

    % z = profile detrended and tapered
    z = profile_taper;
    N = length(z);
    dx = prof_space_m;
    y= fft(z);
    
    % power
    p = y.*conj(y)./(N*dx); 
    
    % put back in dx
    p = p.*dx*dx;
    f = (0:N-1)'/(dx*N);
    p = p(3:N/2);
    f = f(3:N/2);
    
% linear fit to PSD in log space
    if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
        x0_fft = log10(1./f)';
    elseif strcmpi(PSD_xaxis,'frequency') == 1; 
        x0_fft = log10(f)';
    end
        y0_fft = log10(p);
        
     x0_fft = x0_fft(2:end);
     y0_fft = y0_fft(2:end);

    ROLL1 = 1000000000;
    ROLL2 = 0.002304;
    [H_fft,D_fft,RSq_fft,slope_fft,x_fft,yCalc_fft]...
        = fit_lin2log(x0_fft,y0_fft,ROLL1);
    
    H_fft
    D_fft
    RSq_fft;
    
% FIGURE 5 - PSD FFT, EMILY
    xt = 10^4.2;
    yt = 10^-1;
figure(5)
    clf
    if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
        loglog(1./f,p);
        hold on
        xlabel('wavelength (m)')
        set(gca,'fontsize',12,'fontweight','bold')
    elseif strcmpi(PSD_xaxis,'frequency') == 1;
        loglog(f,p);
        hold on
        xlabel('wavenumber (m^{-1})')
        set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')
    end
        loglog(10.^x_fft,10.^yCalc_fft,'r-','linewidth',2)
        
        text_1 = ['Roughness = ',num2str(D_fft)];
        text_2 = ['Hurst = ',num2str(H_fft)];
        text_3 = ['R-squared = ',num2str(RSq_fft)];
        text(xt,yt,{text_1,text_2,text_3},'fontsize',12,'fontweight','bold')
        
        ylabel('power spectral density (m^3)')
        title('power spectral density - FFT') 
        set(gca,'fontsize',12,'fontweight','bold')
        legend('PSD','best fit')
        grid on
       

%% WELCH PSD
 % find the step size of X
    [Pyy,Fw] = pwelch(cum_Disp,[],[],[],fs,'onesided','PSD'); % compute the PSD
  
% linear fit to PSD in log space
   if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
        x0_welch = log10(1./Fw)';
    elseif strcmpi(PSD_xaxis,'frequency') == 1; 
        x0_welch = log10(Fw)';
    end
        y0_welch = log10(Pyy);

    ROLL1 = 1000000000;
    ROLL2 = 0.002304;
    [H_welch,D_welch,RSq_welch,slope_welch,x_welch,yCalc_welch]...
        = fit_lin2log(x0_welch,y0_welch,ROLL1);
    H_welch
    D_welch
    RSq_welch;
    
% FIGURE 4 - PSD WELCH
%     figure(4)
%     plot(log10(F),log10(Pyy)), xlabel('log [\mum^{-1}]'), ylabel('log[PSD[\mum^2/\mum^{-1}]]')
%     title ('Welch PSD')
%     axis tight, grid on
    
    xt = 10^4.2;
    yt = 10^0;
figure(4)
    clf
    if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
        loglog(1./Fw,Pyy);
        hold on
        xlabel('wavelength (m)')
        set(gca,'fontsize',12,'fontweight','bold')
    elseif strcmpi(PSD_xaxis,'frequency') == 1;
        loglog(Fw,Pyy);
        hold on
        xlabel('wavenumber (m^{-1})')
        set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')        
    end
        loglog(10.^x_welch,10.^yCalc_welch,'r-','linewidth',2)
        
        text_1 = ['Roughness = ',num2str(D_welch)];
        text_2 = ['Hurst = ',num2str(H_welch)];
        text_3 = ['R-squared = ',num2str(RSq_welch)];
        text(xt,yt,{text_1,text_2,text_3},'fontsize',12,'fontweight','bold')
        
        ylabel('power spectral density (m^3)')
        title('power spectral density - Welch') 
        legend('PSD','best fit')
        grid on

%% RMS ROUGHNESS
RMS_dev_raw = sqrt(1/length(cum_Disp)*sum(cum_Disp.^2))
sigma_RMS_dt = sqrt(1/length(profile)*sum(profile.^2))
RMS_roughness = rms(cum_Disp)

% Divide fault into windows of varying sizes and calculate RMS error

        figure(9)
            clf
            legend('')
            
    if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        win_size = [5:5:(floor(distalong_transect_m(end)/(1000)))]; % window sizes for sliding window
    else
        win_size = [0.5:0.5:(floor(distalong_transect_m(end)/(1000)))];
    end
    RMS_win = zeros(ones(1,length(win_size)));
    RMS_win_sigma = zeros(ones(1,length(win_size)));
    D_fft_win = zeros(ones(1,length(win_size)));
    H_fft_win = zeros(ones(1,length(win_size)));
    linecolors = jet(length(win_size));
    for jj = 1:length(win_size)
        window = win_size(jj);
        win_end = [window:slide_size:distalong_transect_m(end)/1000];
        win_end = [win_end distalong_transect_m(end)/1000];
        win_start = [0:slide_size:distalong_transect_m(end)/1000];
        num_win = length(win_end);
        win_start = win_start(1:num_win);
        win_center = (win_start+win_end)/2;
        
%         x0_hurst_sort = sort(x0_hurst);
%         rm_infinite = isfinite(x0_hurst_sort);
%         x0_hurst_sort = x0_hurst_sort(rm_infinite);
%         win_start = [x0_hurst_sort(1):slide_size*slide_overlap:x0_hurst_sort(end)];
%         win_end = [win_start(2:end) x0_hurst_sort(end)];
%         num_win = length(win_start);
%         for kk = 1:num_win
%             find_win = (x0_hurst >= win_start(kk) & x0_hurst <= win_end(kk));
%             x0_hurst_win = 
%             y0_hurst_win
%             [Hurst_Thom_slide,Roughness_Thom_slide,RSq_Thom_slide,slope_Thom_slide,...
%                 x_Thom_slide,yCalc_Thom_slide] = fit_lin2log(x0_hurst_win,y0_hurst_win,ROLL1);
%         
        
        RMS_dev_sliding_raw = zeros(ones(1,num_win));
        RMS_dev_sliding_dt = zeros(ones(1,num_win));
        H_fourier_win = zeros(ones(1,num_win));
        D_fourier_win = zeros(ones(1,num_win));
        for ii = 1:num_win
            indeces = ((distalong_transect_m/1000) >= win_start(ii) & (distalong_transect_m/1000) < win_end(ii));
            disp_window = cum_Disp(indeces);

            % raw sigma RMS
            RMS_dev_raw = sqrt(1/length(disp_window)*sum(disp_window.^2));
            RMS_dev_sliding_raw(ii) = RMS_dev_raw;

            % detrended
            meanadj_window = disp_window - mean(disp_window); 
            dt_window = detrend(meanadj_window);
            RMS_dev_detrended = sqrt(1/length(dt_window)*sum(dt_window.^2));
            RMS_dev_sliding_dt(ii) = RMS_dev_detrended;
            
            % FFT - THIBAULT CODE
                %Detrending step:
                disp_window = disp_window - mean(disp_window); %%
                disp_window = detrend(disp_window);
                %Tapering step:
                percent = 0.03;
                disp_window_taper = Ftapering_copy(disp_window,percent);
            % EMILY 
                % z = profile detrended and tapered
                z_win = disp_window_taper;
                N = length(z_win);
                dx = prof_space_m;
                y = fft(z_win);
                p = y.*conj(y)./(N*dx); % power
                % put back in dx
                p = p.*dx*dx;
                f = (0:N-1)'/(dx*N);
                p = p(3:N/2);
                f = f(3:N/2);
            % linear fit to PSD in log space
                if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
                    x0_fft_win = log10(1./f)';
                elseif strcmpi(PSD_xaxis,'frequency') == 1; 
                    x0_fft_win = log10(f)';
                end
                    y0_fft_win = log10(p);

                ROLL1_win = 1000000000;
                [Hurst_fft_win,Roughness_fft_win,RSq_fft_win,slope_fft_win,x_fft_win,yCalc_fft_win]...
                    = fit_lin2log(x0_fft_win,y0_fft_win,ROLL1_win);
                H_fourier_win(ii) = Hurst_fft_win;
                D_fourier_win(ii) = Roughness_fft_win;
                
%                 figure(9)
%                     plot([win_start(ii) win_end(ii)],[sigma_RMS_detrended sigma_RMS_detrended],'-','linewidth',2,'color',linecolors(jj,:))
%                     hold on
%                     legend(strcat('n =',num2str(jj)))
%         
        end
        % sigma_RMS_sliding_dt
        RMS_win(jj) = mean(RMS_dev_sliding_dt);
        H_fft_win(jj) = mean(H_fourier_win);
        D_fft_win(jj) = mean(D_fourier_win);
        RMS_win_sigma(jj) = std(RMS_dev_sliding_dt);  
        
        figure(9)
            plot(win_center,RMS_dev_sliding_dt,'-','linewidth',2,'color',linecolors(jj,:))
            hold on
                    
    end
    figure(9)
        plot(distalong_transect_m(end)/1000,sigma_RMS_dt,'*','MarkerSize',12,'color',linecolors(jj,:))
            tit_text = ['RMS as a function of length scale (window length) for ',EQ_name];
            title(tit_text)
           
        if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
            h(1) = plot(NaN,NaN,'-','color',linecolors(1,:),'linewidth',2);
            h(2) = plot(NaN,NaN,'-','color',linecolors(2,:),'linewidth',2);
            h(3) = plot(NaN,NaN,'-','color',linecolors(5,:),'linewidth',2);
            h(4) = plot(NaN,NaN,'-','color',linecolors(10,:),'linewidth',2);
            h(5) = plot(NaN,NaN,'-','color',linecolors(15,:),'linewidth',2);
            h(6) = plot(NaN,NaN,'-','color',linecolors(20,:),'linewidth',2);
            h(7) = plot(NaN,NaN,'-','color',linecolors(25,:),'linewidth',2);
            h(8) = plot(NaN,NaN,'-','color',linecolors(30,:),'linewidth',2);
            h(9) = plot(NaN,NaN,'-','color',linecolors(35,:),'linewidth',2);
            h(10) = plot(NaN,NaN,'*','color',linecolors(39,:));  
            legend(h,'5 km','10 km','25 km','50 km','75 km','100 km','125 km','150 km','175 km',...
                'full measured fault length');
        elseif strcmpi(EQ,'Kumamoto_LLajoie_40m_XZ.mat') == 1
            h(1) = plot(NaN,NaN,'-','color',linecolors(1,:),'linewidth',2);
            h(2) = plot(NaN,NaN,'-','color',linecolors(2,:),'linewidth',2);
            h(3) = plot(NaN,NaN,'-','color',linecolors(4,:),'linewidth',2);
            h(4) = plot(NaN,NaN,'-','color',linecolors(6,:),'linewidth',2);
            h(5) = plot(NaN,NaN,'-','color',linecolors(8,:),'linewidth',2);
            h(6) = plot(NaN,NaN,'-','color',linecolors(10,:),'linewidth',2);
            h(7) = plot(NaN,NaN,'-','color',linecolors(12,:),'linewidth',2);
            h(8) = plot(NaN,NaN,'-','color',linecolors(14,:),'linewidth',2);
            h(9) = plot(NaN,NaN,'-','color',linecolors(16,:),'linewidth',2);
            h(10) = plot(NaN,NaN,'-','color',linecolors(18,:),'linewidth',2);
            h(11) = plot(NaN,NaN,'-','color',linecolors(20,:),'linewidth',2);
            h(12) = plot(NaN,NaN,'-','color',linecolors(22,:),'linewidth',2);
            h(13) = plot(NaN,NaN,'-','color',linecolors(24,:),'linewidth',2);
            h(14) = plot(NaN,NaN,'-','color',linecolors(26,:),'linewidth',2);
            h(15) = plot(NaN,NaN,'-','color',linecolors(28,:),'linewidth',2);
            h(16) = plot(NaN,NaN,'-','color',linecolors(30,:),'linewidth',2);
            h(17) = plot(NaN,NaN,'*','color',linecolors(30,:));  
            legend(h,'.5 km','1 km','2 km','3 km','4 km','5 km','6 km','7 km','8 km',...
                '9 km','10 km','11 km','12 km','13 km','14 km','15 km','full measured fault length');
        end
        
            xlabel('distance along fault (km)')
            ylabel('RMS deviation of windows (m)')
            set(gca,'fontsize',16,'fontweight','bold') %,'yscale','log','xscale','log')
%             xlim([0 (distalong_transect_m(end)/1000+1)])
%             ylim([0.1 .9])
%             set(gca, 'XTick',[1:15]); % [1 2 3 4 5 6 7 8 9 10 12 15]); % [1:15]); % 
%             set(gca, 'YTick',[0:.1:1])
            grid on
            hold off
    
%% CALCULATE HURST FROM SLOPE OF LOGLOG(wavevector,rms roughness)

if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
    limit1 = 10000000;
    limit2 = [15 100];
    limit3 = [100 160];
        [W,W,W,slope_allpoints,x_RMS_allpoints,y_RMS_fit_allpoints] = fit_lin2log(log10(win_size),log10(RMS_win),limit1);
        [W,W,W,slope_limited,x_RMS_limited,y_RMS_fit_limited] = fit_lin2log(log10(win_size),log10(RMS_win),limit2);
        [W,W,W,slope_limited2,x_RMS_limited2,y_RMS_fit_limited2] = fit_lin2log(log10(win_size),log10(RMS_win),limit3);

        H_RMSslope_allpoints = slope_allpoints
        H_RMSslope_limited = slope_limited
        H_RMSslope_limited2 = slope_limited2
        D_RMSslope_allpoints = 2- slope_allpoints
        D_RMSslope_limited = 2 - slope_limited
        D_RMSslope_limited2 = 2 - slope_limited2
        
elseif strcmpi(EQ,'Kumamoto_LLajoie_40m_XZ.mat') == 1
    limit1 = 10000000;
    limit2 = [.5 10]; 
        [W,W,W,slope_allpoints,x_RMS_allpoints,y_RMS_fit_allpoints] = fit_lin2log(log10(win_size),log10(RMS_win),limit1);
        [W,W,W,slope_limited,x_RMS_limited,y_RMS_fit_limited] = fit_lin2log(log10(win_size),log10(RMS_win),limit2);
        H_RMSslope_allpoints = slope_allpoints
        H_RMSslope_limited = slope_limited
        D_RMSslope_allpoints = 2- slope_allpoints
        D_RMSslope_limited = 2 - slope_limited
elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
    limit1 = 10000000;
    limit2 = [1.5 50]; 
        [W,W,W,slope_allpoints,x_RMS_allpoints,y_RMS_fit_allpoints] = fit_lin2log(log10(win_size),log10(RMS_win),limit1);
        [W,W,W,slope_limited,x_RMS_limited,y_RMS_fit_limited] = fit_lin2log(log10(win_size),log10(RMS_win),limit2);
        H_RMSslope_allpoints = slope_allpoints
        H_RMSslope_limited = slope_limited
        D_RMSslope_allpoints = 2- slope_allpoints
        D_RMSslope_limited = 2 - slope_limited
elseif strcmpi(EQ,'HectorMine_Milliner.csv') == 1
    limit1 = 10000000;
    limit2 = [0 14];
    limit3 = [14 36];
        [W,W,W,slope_allpoints,x_RMS_allpoints,y_RMS_fit_allpoints] = fit_lin2log(log10(win_size),log10(RMS_win),limit1);
        [W,W,W,slope_limited,x_RMS_limited,y_RMS_fit_limited] = fit_lin2log(log10(win_size),log10(RMS_win),limit2);
        [W,W,W,slope_limited2,x_RMS_limited2,y_RMS_fit_limited2] = fit_lin2log(log10(win_size),log10(RMS_win),limit3);

        H_RMSslope_allpoints = slope_allpoints
        H_RMSslope_limited = slope_limited
        H_RMSslope_limited2 = slope_limited2
        D_RMSslope_allpoints = 2- slope_allpoints
        D_RMSslope_limited = 2 - slope_limited
        D_RMSslope_limited2 = 2 - slope_limited2
end

%% SAVE DATA
    if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        save('Balochistan_plot_data.mat','win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','x_RMS_limited2','y_RMS_fit_limited2',...
            'H_RMSslope_allpoints','H_RMSslope_limited','H_RMSslope_limited2')
    elseif strcmpi(EQ,'Kumamoto_LLajoie_40m_XZ.mat') == 1
        save('Kumamoto_plot_data.mat','win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','H_RMSslope_allpoints',...
            'H_RMSslope_limited')
    elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
        save('Landers_plot_data.mat','win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','H_RMSslope_allpoints',...
            'H_RMSslope_limited')
    elseif strcmpi(EQ,'HectorMine_Milliner.csv') == 1
        save('HectorMine_plot_data.mat','win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','x_RMS_limited2','y_RMS_fit_limited2',...
            'H_RMSslope_allpoints','H_RMSslope_limited','H_RMSslope_limited2')
    end
        
%% RMS,D, H v length - log
    figure(6)
    clf
    plot(win_size,H_fft_win,'b-','linewidth',2)
    hold on
    plot(distalong_transect_m(end)/1000,H_fft,'b*','MarkerSize',12)
    hold on
    plot(win_size,D_fft_win,'r-','linewidth',2)
    hold on
    plot(distalong_transect_m(end)/1000,D_fft,'r*','MarkerSize',12)
    hold on
    plot(win_size,RMS_win,'g-','linewidth',2)
    hold on
    plot(distalong_transect_m(end)/1000,sigma_RMS_dt,'g*','MarkerSize',12)
    hold on
    errorbar(win_size,RMS_win,RMS_win_sigma,'color','g')
    
    if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        xlim([4 200])
        ylim([0.02 3.1])
%         set(gca, 'XTick',[0:15]); % 0 1 2 3 4 5 6 7 8 9 10 12 15])
%         set(gca, 'YTick',[0:.1:1])
    elseif strcmpi(EQ,'Kumamoto_LLajoie_40m_XZ.mat') == 1
        hold on
        plot(40/1000,0.23,'go')
        xlim([40/1000 (distalong_transect_m(end)/1000+1)])
        ylim([-1 2.5])
%         set(gca, 'XTick',[0:15]); % 0 1 2 3 4 5 6 7 8 9 10 12 15])
%         set(gca, 'YTick',[0:.1:1])
    elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
        xlim([.4 (distalong_transect_m(end)/1000+1)])
        ylim([.02 2])
    end

        tit_text = ['RMS, D, H as a function of length scale (window length) for ',EQ_name];
        title(tit_text)
        xlabel('length of window for calculation (km)')
        ylabel('average value of roughness parameter (RMS in m)')
        legend('Hurst (H) for windowed data','Hurst (H) for full dataset','Fractal dim (D) for windowed data',...
            'Fractal dim (D) for full dataset','RMS for windowed data','RMS for full dataset',...
            '1-sigma error bars for RMS','dataset noise level')
        set(gca,'fontsize',16,'fontweight','bold','yscale','log','xscale','log')
        grid on
        
%% RMS,D, H v length - linspace
    figure(7)
    clf
    plot(win_size,H_fft_win,'b-','linewidth',2)
    hold on
    plot(distalong_transect_m(end)/1000,H_fft,'b*','MarkerSize',12)
    hold on
    plot(win_size,D_fft_win,'r-','linewidth',2)
    hold on
    plot(distalong_transect_m(end)/1000,D_fft,'r*','MarkerSize',12)
    hold on
    plot(win_size,RMS_win,'g-','linewidth',2)
    hold on
    plot(distalong_transect_m(end)/1000,sigma_RMS_dt,'g*','MarkerSize',12)
    hold on
    errorbar(win_size,RMS_win,RMS_win_sigma,'color','g')
        tit_text = ['H and D as a function of length scale for ',EQ_name];
        title(tit_text)
        xlabel('length of window for calculation (km)')
        ylabel('average value of roughness parameter')
        legend('Hurst (H) for windowed data','Hurst (H) for full dataset','Fractal dim (D) for windowed data',...
            'Fractal dim (D) for full dataset','RMS for windowed data','RMS for full dataset',...
            '1-sigma error bars for RMS')
        set(gca,'fontsize',16,'fontweight','bold')
        grid on
        
%% FIGURE 8 - HURST FROM RMS ROUGHNESS
    figure(8)
    clf
    plot(win_size,RMS_win,'g-','linewidth',2)
    hold on
    plot(distalong_transect_m(end)/1000,sigma_RMS_dt,'g*','MarkerSize',12)
    hold on
    errorbar(win_size,RMS_win,RMS_win_sigma,'color','g')
    hold on
%     plot(10.^(x_RMS_allpoints),10.^(y_RMS_fit_allpoints),'k--','linewidth',2)
%     hold on
    plot(10.^(x_RMS_limited),10.^(y_RMS_fit_limited),'k:','linewidth',2)
    
    if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        hold on
        plot(10.^(x_RMS_limited2),10.^(y_RMS_fit_limited2),'k-.','linewidth',2)
        xlim([4 200])
        ylim([0.08 3.1])
%         set(gca, 'XTick',[0:15]); % 0 1 2 3 4 5 6 7 8 9 10 12 15])
%         set(gca, 'YTick',[0:.1:1])
    elseif strcmpi(EQ,'Kumamoto_LLajoie_40m_XZ.mat') == 1
        hold on
        plot(40/1000,0.23,'go')
        xlim([40/1000 (distalong_transect_m(end)/1000+1)])
        ylim([-1.5 10^-.1])
%         set(gca, 'XTick',[0:15]); % 0 1 2 3 4 5 6 7 8 9 10 12 15])
%         set(gca, 'YTick',[0:.1:1])
    elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
        xlim([.4 (distalong_transect_m(end)/1000+1)])
        ylim([.1 2])
    elseif strcmpi(EQ,'HectorMine_Milliner.csv') == 1
        hold on
        plot(10.^(x_RMS_limited2),10.^(y_RMS_fit_limited2),'k-.','linewidth',2)
    end

        tit_text = ['RMS as a function of window size for ',EQ_name];
        title(tit_text)
        xlabel('length of window for calculation (km)')
        ylabel('average value of roughness parameter (RMS in m)')
        tstring1 = ['linear fit to RMS (all points): slope = H = ',num2str(H_RMSslope_allpoints)];
        tstring2 = ['linear fit to RMS (limited range): slope = H = ',num2str(H_RMSslope_limited)];
        if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
            tstring3 = ['linear fit to RMS (limited upper range): slope = H = ',num2str(H_RMSslope_limited2)];
            legend('RMS for windowed data','RMS for full dataset',...
                '1-sigma error bars for RMS',tstring2,tstring3)
        elseif strcmpi(EQ,'HectorMine_Milliner.csv') == 1
            tstring3 = ['linear fit to RMS (limited upper range): slope = H = ',num2str(H_RMSslope_limited2)];
            legend('RMS for windowed data','RMS for full dataset',...
                '1-sigma error bars for RMS',tstring2,tstring3)
        else
            legend('RMS for windowed data','RMS for full dataset',...
                '1-sigma error bars for RMS',tstring2)
        end
        set(gca,'fontsize',16,'fontweight','bold','yscale','log','xscale','log')
        grid on
        

%% TEST PLOT
    figure(1)
    clf
    plot(distalong_transect_km,disp_m,'b-','linewidth',2)
    hold on
    plot(distalong_transect_km,profile,'g-','linewidth',2)
    hold on
    plot(distalong_transect_km,profile_taper,'r-','linewidth',2)
    hold on
%     for ii = 1:div_fault
%         x_lim = [distalong_transect_km(seg_ind_start(ii)) distalong_transect_km(seg_ind_end(ii))];
%         plot(x_lim,[sigma_RMS_sliding_raw(ii) sigma_RMS_sliding_raw(ii)],'b--')
%         hold on
%         plot(x_lim,[sigma_RMS_sliding_dt(ii) sigma_RMS_sliding_dt(ii)],'b--','linewidth',2)
%         hold on
%     end
        tit_text = ['Slip distributions for ',EQ_name];
        title(tit_text)
        xlabel('distance along profile (km)')
        ylabel('lateral offset (m)')
        legend('raw slip','detrended slip','detrended and tapered slip')
        set(gca,'fontsize',16,'fontweight','bold')
        grid on

% % figure 25
%     Figures25_26_ProfileStatsResults_Kum