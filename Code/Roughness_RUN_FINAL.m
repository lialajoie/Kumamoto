% Lia Lajoie
% CSM

% Calculate rouchnesses and make plots
clear all

%% USER INPUTS
PSD_xaxis = 'wavelength'; % 'frequency' or 'wavelength'
slide_size = .5; % in km - how much windows slide
define_win = [3 5 10]; % sliding window sizes for strain and strain gradient
strain_plot = 1; %1, 2, 3, 4 correspond to sliding windows 1, 3, 5, 10, respectively

KUM_FAULT = 8; % for Kumamoto only! Choose 1-6 individual fault segments, to run individually, or anter any other name/number into here to run full-fault
 
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
   % EQ = 'Kumamoto_40m_3seg_AFTERBUG.mat';
   % EQ = 'Landers_disp_data_no_header.csv'; % from Milliner
   % EQ = 'HectorMine_Milliner.csv'; % from Milliner
   % EQ = 'Balochistan_Vallage_FarField.txt'; % from Vallage
   % EQ = 'Palu_long_disp.txt'; % Palu from Socquet
   EQ = 'Baluchistan_stack_out_nutm41.txt'; % from Zinke
   % EQ = 'EMC1';
   % EQ = 'EMC2';
   % EQ = 'EMC3';
   % EQ = 'EMC4';
   %  EQ = 'EMC5';
    % EQ = 'EMC8';
    
%% IMPORT KUMAMOTO FIELD DATA FOR GENERALIZED 3 SEGMENT FAULT
load('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/FieldData/KumFieldData_plot_FINAL.mat')
f_dist_km = Field.F_distalongfault/1000;
f_Rslip = Field.f_Rslip;

%% RETRIEVE EARTHQUAKE DATA
    EQ_dir = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/';
    EQ_file = [EQ_dir,EQ];
    
    if strcmpi(EQ,'Kumamoto_XZ.mat') == 1
        S = open(EQ_file);
        distalong_transect_m = S.x;
        cumdisp_m = S.z;
        EQ_name = 'Kumamoto - Candela';
    elseif strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
        load(EQ_file);
        prof_space = 40;
        
        trans_bend1 = 2.7597; % km along profile
        trans_bend2 = 9.438 + trans_bend1; %km along profile
        
        % flip all dist along transect LR
        prof_x = cell2mat({ProfsALL.ctr_easting});
        prof_y = cell2mat({ProfsALL.ctr_northing});  
        cumdisp_m = cell2mat({ProfsALL.linear_disp_cum});  
        err = cell2mat({ProfsALL.linear_err_cum});
        
        EQ_name = 'Kumamoto - Lajoie';
        distalong_transect_m = fliplr([0:length(prof_x)-1])*prof_space;
        distalong_transect_km = distalong_transect_m/1000;
        
            % individual scarp data
                Scarps_x = [];
                Scarps_y = [];
                Scarps_distalongtransect = [];
                Scarps_lat = [];
                Scarps_lat_err = [];
                Scarps_cumLat = [];
                Scarps_cumLat_err = [];
                Cumdisp_m = zeros(size(prof_x));
                for ii = 1:length(prof_x)
                        Scarps_alongprof = cell2mat({ProfsALL(ii).scarp_alongprof});
                        % calculate profile parameters
                        P_endx = cell2mat({ProfsALL(ii).endpts_x});
                        P_endy = cell2mat({ProfsALL(ii).endpts_y});
                        dx = (P_endx(2)-P_endx(1));
                        dy = (P_endy(2)-P_endy(1));
                        prof_len = sqrt(dx^2+dy^2)/2;
                        
                        Scarp_vals = cell2mat({ProfsALL(ii).linear_extrap_disp});
                        sum_scarps = sum(Scarp_vals);
                        d_alongtran = ones(size(Scarp_vals))'*distalong_transect_km(ii);
                        Cumdisp_m(ii) = sum_scarps;

                        scarps_x = zeros(1,length(Scarps_alongprof));
                        scarps_y = scarps_x;
                        for jj = 1:length(Scarps_alongprof)
                            scarp_dist = Scarps_alongprof(jj);
                            a = scarp_dist/(prof_len*2);
                            scarp_x = prof_x(ii)+dx*a;
                            scarp_y = prof_y(ii)+dy*a;
                            scarps_x(jj) = scarp_x;
                            scarps_y(jj) = scarp_y;

                            dist_scarp_tocen = sqrt((prof_x(ii)-scarp_x)^2+(prof_y(ii)-scarp_y)^2);
                        end
                        scarps_lat = Scarp_vals;
                        scarps_lat_err = cell2mat({ProfsALL(ii).linear_extrap_error});
                        scarps_cumLat = cell2mat({ProfsALL(ii).linear_disp_cum});
                        scarps_cumLat_err = cell2mat({ProfsALL(ii).linear_err_cum});

                %         figure(3)
                %             clf
                %             plot(P_endx,P_endy,'k:')
                %             hold on
                %             scatter(Profs_ctr_x(ii),Profs_ctr_y(ii),10,'g','filled')
                %             hold on
                %             scatter(scarps_x,scarps_y,5,'r','filled')

                        % save in vector
                        Scarps_x = [Scarps_x scarps_x];
                        Scarps_y = [Scarps_y scarps_y];
                        Scarps_distalongtransect = [Scarps_distalongtransect d_alongtran];
                        Scarps_lat = [Scarps_lat scarps_lat'];
                        Scarps_lat_err = [Scarps_lat_err scarps_lat_err'];
                        Scarps_cumLat = [Scarps_cumLat scarps_cumLat'];
                        Scarps_cumLat_err = [Scarps_cumLat_err scarps_cumLat_err'];
                end
    
        % pseudo- full-slip
%         add_dist = distalong_transect_m + max(distalong_transect_m) + (distalong_transect_m(2)-distalong_transect_m(1));
%         distalong_transect_m = [distalong_transect_m add_dist];
%         disp_m = [disp_m fliplr(disp_m)
    elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
        S = importdata(EQ_file);
        distalong_transect_km = S(:,1);
        distalong_transect_m = distalong_transect_km*1000;
        cumdisp_m = S(:,2);
        EQ_name = 'Landers - Milliner';
    elseif strcmpi(EQ,'HectorMine_Milliner.csv') == 1
        S = importdata(EQ_file);
        cumdisp_m = S(:,1);
        distalong_transect_m = [0:1:length(cumdisp_m)-1]*100;
        EQ_name = 'Hector Mine - Milliner';
    elseif strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        S = importdata(EQ_file);
        distalong_transect_km = S(:,1); % distance from south
        distalong_transect_km = sort(distalong_transect_km);
        distalong_transect_m = distalong_transect_km*1000;
        cumdisp_m = S(:,2); % fault parallel
        EQ_name = 'Balochistan - Vallage';
%     elseif isempty(regexp(EQ,'EMC')) ~= 1
%         EQ_dir = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/DipResults_editStrike_Fault';
%         EQ_file = [EQ_dir,EQ(4),'.mat'];
%         S = open(EQ_file);
%         disp_m = S.mag_disp;
%         distalong_transect_m = [0:length(disp_m)-1]*300;
%         EQ_name = EQ;
    elseif strcmpi(EQ,'Baluchistan_stack_out_nutm41.txt') == 1
        S = importdata(EQ_file);
        cumdisp_m = S(:,3);
        UTM_east = S(:,1);
        UTM_north = S(:,2);
        dist_m = zeros(size(UTM_east));
        for ii = 1:length(UTM_east)-1
            dist_m(ii+1) = sqrt((UTM_east(ii+1)- UTM_east(ii))^2 + (UTM_north(ii+1)- UTM_north(ii))^2);
        end
%         dist_m = dist_km*1000;
        distalong_transect_m = cumsum(dist_m);
        EQ_name = 'Balochistan - Zinke';
    elseif strcmpi(EQ,'Palu_long_disp.txt') == 1
        S = importdata(EQ_file);
        cumdisp_m = S(:,5);
        UTM_east = S(:,3);
        UTM_north = S(:,4);
        dist_m = zeros(size(UTM_east));
        for ii = 1:length(UTM_east)-1
            dist_m(ii+1) = sqrt((UTM_east(ii+1)- UTM_east(ii))^2 + (UTM_north(ii+1)- UTM_north(ii))^2);
        end
        distalong_transect_m = cumsum(dist_m);
        EQ_name = 'Palu - Socquet';
    elseif strcmpi(EQ,'EMC5') == 1
        EQ_file4 = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/DipResults_editStrike_Fault4.mat';
        S4 = open(EQ_file4);
        disp_m4 = S4.mag_disp;
        EQ_file5 = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/DipResults_editStrike_Fault5.mat';
        S5 = open(EQ_file5);
        disp_m5 = S5.mag_disp;
        
        cumdisp_m = [disp_m5 disp_m4];
        distalong_transect_m = [0:length(cumdisp_m)-1]*300;
        EQ_name = 'EMC Paso Superior fault';
    else
        S = open(EQ_file);
        distalong_transect_m = S.X;
        cumdisp_m = S.Z;
        EQ_name = 'Tibo EQ';
    end
    distalong_transect_km = distalong_transect_m/1000;
    
%% LOAD FAULT DATA
if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
    load('Fault_segments.mat')
    Faults(1).xy = Fault1; Faults(2).xy = Fault2; Faults(3).xy = Fault3;
    Faults(4).xy = Fault4; Faults(5).xy = Fault5; Faults(6).xy = Fault6;

        for mm = 1:max(size(Faults))
            ind1 = [];
            for kk = 1:max(size(Faults(mm).xy))
                ind = find(Faults(mm).xy(kk,1) == Scarps_x);
                ind1 = [ind1 ind];
            end
            Faults(mm).ind = ind1;
        end

            Faults(1).color = [102 51 153]/256;
            Faults(2).color = [255 102 0]/256;
            Faults(3).color = [0 0 1];
            Faults(4).color = [51 153 0]/256;
            Faults(5).color = [1, 0, 1];
            Faults(6).color = [0 153 204]/256;
end

%% CALCULATE PARAMS
prof_space_m = abs(mean(distalong_transect_m(2:end)-distalong_transect_m(1:end-1)));
std_profspace = std(distalong_transect_m(2:end)-distalong_transect_m(1:end-1));
Si = prof_space_m; %/1000; %/1000; % sampling interval  = profile spacing in km 
xmin = (prof_space_m*2);
cum_Disp = cumdisp_m;
N = length(cum_Disp);

%% STRAIN AND STRAIN GRADIENT - WHOLE FAULT
% strain = ds/dx; strain grad = d2s/dx2
    
% point to point strain 
ds = cum_Disp(2:end) - cum_Disp(1:end-1);
dx = distalong_transect_m(2:end)-distalong_transect_m(1:end-1);
x_strain = (distalong_transect_m(2:end)+distalong_transect_m(1:end-1))/2;
strain = ds./dx;

% point to point strain gradient
d2s = strain(2:end) - strain(1:end-1);
dx2 = x_strain(2:end)-x_strain(1:end-1);
x_sgrad = (x_strain(2:end)+x_strain(1:end-1))/2;
strain_grad = d2s./dx2;

% store data
Slide(1).strain = strain;
Slide(1).strain_std = zeros(size(strain));
Slide(1).x_strain = x_strain;
Slide(1).strain_grad = strain_grad;
Slide(1).strain_grad_std = zeros(size(strain_grad));
Slide(1).x_sgrad = x_sgrad;
       
% sliding average of strain and strain gradient
for mm = 1:length(define_win);
    slide_win = define_win(mm);
    
    % SLIDING WINDOWS STRAIN
    slide_str_end = [slide_win:length(strain)];
    slide_str_start = [1:length(strain)];
    slide_str_num = length(slide_str_end);
    slide_str_start = slide_str_start(1:slide_str_num);
    slide_str_center = (distalong_transect_m(slide_str_start) + distalong_transect_m(slide_str_end))/2;
    
    strain_slide = zeros(ones(1,slide_str_num));
    strain_std_slide = zeros(ones(1,slide_str_num));
    for nn = 1:slide_str_num
        str = strain(slide_str_start(nn):slide_str_end(nn));
        mean_str = mean(str);
        std_str = std(str);
        
        strain_slide(nn) = mean_str;
        strain_std_slide(nn) = std_str;
    end
    
    % save data
    Slide(mm + 1).strain = strain_slide;
    Slide(mm + 1).strain_std= strain_std_slide;
    Slide(mm + 1).x_strain = slide_str_center;
    
    % SLIDING WINDOWS STRAIN GRADIENT
    slide_grad_end = [slide_win:length(strain_grad)];
    slide_grad_start = [1:length(strain_grad)];
    slide_grad_num = length(slide_grad_end);
    slide_grad_start = slide_grad_start(1:slide_grad_num);
    slide_grad_center = (x_sgrad(slide_grad_start) + x_sgrad(slide_grad_end))/2;

    sgrad_slide = zeros(ones(1,slide_grad_num));
    sgrad_std_slide = zeros(ones(1,slide_grad_num));
    for nn = 1:slide_grad_num
        sgrad = strain_grad(slide_grad_start(nn):slide_grad_end(nn));
        mean_sgrad = mean(sgrad);
        std_sgrad = std(sgrad);
        
        sgrad_slide(nn) = mean_sgrad;
        sgrad_std_slide(nn) = std_sgrad;
    end
    
    % save data
    Slide(mm + 1).strain_grad = sgrad_slide;
    Slide(mm + 1).strain_grad_std= sgrad_std_slide;
    Slide(mm + 1).x_sgrad = slide_grad_center;
end

%% STRAIN AND STRAIN GRADIENT - INDIVIDUAL FAULT SEGMENTS
% strain = ds/dx; strain grad = d2s/dx2

if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
    for yy = 1:max(size(Faults))
        % empty vectors
        x_strain = [];
        x_sgrad = [];

        % point to point strain 
        disp = Scarps_lat(Faults(yy).ind);
        ds = disp(1:end-1) - disp(2:end);
        dx = zeros(size(ds));
        dx_x = dx;
        dx_y = dy;
        for zz1 = 1:length(dx)
            X1 = Faults(yy).xy(zz1,1); X2 = Faults(yy).xy(zz1+1,1);
            Y1 = Faults(yy).xy(zz1,2); Y2 = Faults(yy).xy(zz1+1,2);
            dx(zz1) = sqrt((X1-X2)^2 + (Y1-Y2)^2);
            dx_x(zz1) = (X1+X2)/2;
            dx_y(zz1) = (Y1+Y2)/2;
        end
        x_strain = (Scarps_distalongtransect(Faults(yy).ind(2:end))+Scarps_distalongtransect(Faults(yy).ind(1:end-1)))/2;
        strain = ds./dx;

        % point to point strain gradient
        d2s = strain(1:end-1) - strain(2:end);
        dx2 = zeros(size(d2s));
        for zz2 = 1:length(dx2)
            XX1 = dx_x(zz2); XX2 = dx_x(zz2+1);
            YY1 = dx_y(zz2); YY2 = dx_y(zz2+1);
            dx2(zz2) = sqrt((XX1-XX2)^2 + (YY1-YY2)^2);
        end
        x_sgrad = (x_strain(2:end)+x_strain(1:end-1))/2;
        strain_grad = d2s./dx2;

        % store data
        Faults(yy).Slide(1).strain = strain;
        Faults(yy).Slide(1).dx = dx;
        Faults(yy).Slide(1).strain_std = zeros(size(strain));
        Faults(yy).Slide(1).x_strain = x_strain;
        Faults(yy).Slide(1).strain_grad = strain_grad;
        Faults(yy).Slide(1).strain_grad_std = zeros(size(strain_grad));
        Faults(yy).Slide(1).x_sgrad = x_sgrad;

        % sliding average of strain and strain gradient
        for mm = 1:length(define_win);
            slide_win = define_win(mm);

            % SLIDING WINDOWS STRAIN
            slide_str_end = [slide_win:length(strain)];
            slide_str_start = [1:length(strain)];
            slide_str_num = length(slide_str_end);
            slide_str_start = slide_str_start(1:slide_str_num);
            slide_str_center = (x_strain(slide_str_start) + x_strain(slide_str_end))/2;

            strain_slide = zeros(ones(1,slide_str_num));
            strain_std_slide = zeros(ones(1,slide_str_num));
            for nn = 1:slide_str_num
                str = strain(slide_str_start(nn):slide_str_end(nn));
                mean_str = mean(str);
                std_str = std(str);

                strain_slide(nn) = mean_str;
                strain_std_slide(nn) = std_str;
            end

            % save data
            Faults(yy).Slide(mm + 1).strain = strain_slide;
            Faults(yy).Slide(mm + 1).strain_std= strain_std_slide;
            Faults(yy).Slide(mm + 1).x_strain = slide_str_center;

            % SLIDING WINDOWS STRAIN GRADIENT
            slide_grad_end = [slide_win:length(strain_grad)];
            slide_grad_start = [1:length(strain_grad)];
            slide_grad_num = length(slide_grad_end);
            slide_grad_start = slide_grad_start(1:slide_grad_num);
            slide_grad_center = (x_sgrad(slide_grad_start) + x_sgrad(slide_grad_end))/2;

            sgrad_slide = zeros(ones(1,slide_grad_num));
            sgrad_std_slide = zeros(ones(1,slide_grad_num));
            for nn = 1:slide_grad_num
                sgrad = strain_grad(slide_grad_start(nn):slide_grad_end(nn));
                mean_sgrad = mean(sgrad);
                std_sgrad = std(sgrad);

                sgrad_slide(nn) = mean_sgrad;
                sgrad_std_slide(nn) = std_sgrad;
            end

            % save data
            Faults(yy).Slide(mm + 1).strain_grad = sgrad_slide;
            Faults(yy).Slide(mm + 1).strain_grad_std= sgrad_std_slide;
            Faults(yy).Slide(mm + 1).x_sgrad = slide_grad_center;
        end
    end
end

%% LIMIT FAULT
if KUM_FAULT == 1 || KUM_FAULT == 2 || KUM_FAULT == 3 || ...
        KUM_FAULT == 4 || KUM_FAULT == 5 || KUM_FAULT == 6;
    DISP = Scarps_lat(Faults(KUM_FAULT).ind);
    DIST = Scarps_distalongtransect(Faults(KUM_FAULT).ind);
else
    DISP = cum_Disp;
    DIST = distalong_transect_m/1000;
end

    %% PSD - THOMSON
    x = DISP;
    fs = 1/Si; % sampling frequency
    [pxx,f] = pmtm(x,[],length(x),fs); % multitaper PSD
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
            [H_Thom,D_Thom,RSq_Thom,slope_Thom,x_Thom,yCalc_Thom]...
                = fit_lin2log(x0_hurst,y0_hurst,ROLL1);

        H_Thom
        D_Thom
        RSq_Thom;

    %% PSD - FFT

    % THIBAULT CODE
        %Detrending step:
        profile = DISP;
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

    %% WELCH PSD
     % find the step size of X
        [Pyy,Fw] = pwelch(DISP,[],[],[],fs,'onesided','PSD'); % compute the PSD

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

    %% RMS ROUGHNESS
    RMS_dev_raw = sqrt(1/length(DISP)*sum(DISP.^2))
    sigma_RMS_dt = sqrt(1/length(profile)*sum(profile.^2))
    % RMS_roughness = rms(cum_Disp)
    
    % Divide fault into windows of varying sizes and calculate RMS error

    figure(5)
        clf
        legend('')
        win_length = max(DIST) - min(DIST);

        if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
            win_size = [5:5:floor(max(DIST))]; % window sizes for sliding window
        elseif strcmpi(EQ,'Palu_long_disp.txt') == 1
            win_size = [5:5:floor(max(DIST))]; % window sizes for sliding window
        elseif isempty(regexp(EQ,'EMC')) ~= 1
            win_size = [2:2:floor(max(DIST))];
        elseif strcmpi(EQ,'Baluchistan_stack_out_nutm41.txt') == 1
            win_size = [50:50:floor(max(DIST))]; 
        else
            win_size = [0.5:0.5:floor(win_length)];
        end

        RMS_win = zeros(ones(1,length(win_size)));
        RMS_win_sigma = zeros(ones(1,length(win_size)));
        D_fft_win = zeros(ones(1,length(win_size)));
        H_fft_win = zeros(ones(1,length(win_size)));
        linecolors = jet(length(win_size));
        mean_disp_win = zeros(ones(1,length(win_size)));
        R_RMSdisp = zeros(ones(1,length(win_size)));
        P_RMSdisp = zeros(ones(1,length(win_size)));
        R_StrDisp = zeros(ones(1,length(win_size)));
        P_StrDisp = zeros(ones(1,length(win_size)));
        R_StrRMS = zeros(ones(1,length(win_size)));
        P_StrRMS = zeros(ones(1,length(win_size)));
        lags_max = zeros(ones(1,length(win_size)));
        for jj = 1:length(win_size)
            window = win_size(jj);
            if min(DIST) > min(win_size)
                win_end = [(min(DIST)+window):slide_size:max(DIST)];
                win_start = [min(DIST):slide_size:max(DIST)];
            else
                win_end = [window:slide_size:max(DIST)];
                win_start = [0:slide_size:max(DIST)];
            end
            win_end = [win_end max(DIST)];
            num_win = length(win_end);
            win_start = win_start(1:num_win);
            win_center = (win_start+win_end)/2;

            RMS_dev_sliding_raw = zeros(ones(1,num_win));
            RMS_dev_sliding_dt = zeros(ones(1,num_win));
            H_fourier_win = zeros(ones(1,num_win));
            D_fourier_win = zeros(ones(1,num_win));
            mean_disp = zeros(ones(1,num_win));
            mean_str_win = zeros(ones(1,num_win));
            for ii = 1:num_win
                indeces = ((DIST) >= win_start(ii) & (DIST) < win_end(ii));
                disp_window = DISP(indeces);
                mean_disp(ii) = mean(disp_window);

                % strain
                ds_win = disp_window(2:end) - disp_window(1:end-1);
                str_win = ds_win/40;
                mean_str = mean(abs(str_win));
                mean_str_win(ii) = mean_str;

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
                    p = p(3:ceil(N/2));
                    f = f(3:ceil(N/2));
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

%                     figure(5)
%                         plot([win_start(ii) win_end(ii)],[sigma_RMS_detrended sigma_RMS_detrended],'-','linewidth',2,'color',linecolors(jj,:))
%                         hold on
%                         legend(strcat('n =',num2str(jj)))
    %         
            end
            
            if num_win > 3
                [R_rms,P_rms] = corrcoef(mean_disp,RMS_dev_sliding_dt);
                [r,lags] = xcorr(mean_disp,RMS_dev_sliding_dt);
                max_r = find(r == max(r));
                max_lag = lags(max_r);

                [R_str,P_str] = corrcoef(mean_disp,mean_str_win);
                [R_strRMS,P_strRMS] = corrcoef(RMS_dev_sliding_dt,mean_str_win);
            else 
                R_rms = [NaN NaN]; P_rms= [NaN NaN];
                R_str = [NaN NaN]; P_str = [NaN NaN];
                R_strRMS = [NaN NaN]; P_strRMS = [NaN NaN];
                max_lag = NaN;
            end

            % save data
            RMS_win(jj) = mean(RMS_dev_sliding_dt);
            H_fft_win(jj) = mean(H_fourier_win);
            D_fft_win(jj) = mean(D_fourier_win);
            RMS_win_sigma(jj) = std(RMS_dev_sliding_dt);  
            R_RMSdisp(jj) = R_rms(2);
            P_RMSdisp(jj) = P_rms(2);
            R_StrDisp(jj) = R_str(2);
            P_StrDisp(jj) = P_str(2);
            R_StrRMS(jj) = R_strRMS(2);
            P_StrRMS(jj) = P_strRMS(2);
            lags_max(jj) = max_lag;

            % save data for specific fault lengths
            if window == 1
                RMS_1km = RMS_dev_sliding_dt;
                mean_slip_1km = mean_disp;
                win_center_1km = win_center;
                win_start_1km = win_start;
                win_end_1km = win_end;
                mean_str_1km = mean_str_win;
            elseif window == 3
                RMS_3km = RMS_dev_sliding_dt;
                mean_slip_3km = mean_disp;
                win_center_3km = win_center;
                win_start_3km = win_start;
                win_end_3km = win_end;
                mean_str_3km = mean_str_win;
            elseif window == 5
                RMS_5km = RMS_dev_sliding_dt;
                mean_slip_5km = mean_disp;
                win_center_5km = win_center;
                win_start_5km = win_start;
                win_end_5km = win_end;
                mean_str_5km = mean_str_win;
            end


            figure(5)
                plot(win_center,RMS_dev_sliding_dt,'-','linewidth',2,'color',linecolors(jj,:))
                hold on

        end

        figure(5)
            plot(max(DIST),sigma_RMS_dt,'*','MarkerSize',12,'color',linecolors(jj,:))
            if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1 % plot inflection points in transect
                hold on
                plot([trans_bend1 trans_bend1],[0 .9],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                hold on
                plot([trans_bend2 trans_bend2],[0 .9],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
            end
                tit_text = ['RMS as a function of length scale (window length) for ',EQ_name];
                title(tit_text)

            if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
                h2(1) = plot(NaN,NaN,'-','color',linecolors(1,:),'linewidth',2);
                h2(2) = plot(NaN,NaN,'-','color',linecolors(2,:),'linewidth',2);
                h2(3) = plot(NaN,NaN,'-','color',linecolors(5,:),'linewidth',2);
                h2(4) = plot(NaN,NaN,'-','color',linecolors(10,:),'linewidth',2);
                h2(5) = plot(NaN,NaN,'-','color',linecolors(15,:),'linewidth',2);
                h2(6) = plot(NaN,NaN,'-','color',linecolors(20,:),'linewidth',2);
                h2(7) = plot(NaN,NaN,'-','color',linecolors(25,:),'linewidth',2);
                h2(8) = plot(NaN,NaN,'-','color',linecolors(30,:),'linewidth',2);
                h2(9) = plot(NaN,NaN,'-','color',linecolors(35,:),'linewidth',2);
                h2(10) = plot(NaN,NaN,'*','color',linecolors(39,:));  
                legend(h2,'5 km','10 km','25 km','50 km','75 km','100 km','125 km','150 km','175 km',...
                    'full measured fault length');
            elseif strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
                if KUM_FAULT == 1 || KUM_FAULT == 2 || KUM_FAULT == 3 || ...
                    KUM_FAULT == 4 || KUM_FAULT == 5 || KUM_FAULT == 6;
                
                else
                    h2(1) = plot(NaN,NaN,'-','color',linecolors(1,:),'linewidth',2);
                    h2(2) = plot(NaN,NaN,'-','color',linecolors(2,:),'linewidth',2);
                    h2(3) = plot(NaN,NaN,'-','color',linecolors(4,:),'linewidth',2);
                    h2(4) = plot(NaN,NaN,'-','color',linecolors(6,:),'linewidth',2);
                    h2(5) = plot(NaN,NaN,'-','color',linecolors(8,:),'linewidth',2);
                    h2(6) = plot(NaN,NaN,'-','color',linecolors(10,:),'linewidth',2);
                    h2(7) = plot(NaN,NaN,'-','color',linecolors(12,:),'linewidth',2);
                    h2(8) = plot(NaN,NaN,'-','color',linecolors(14,:),'linewidth',2);
                    h2(9) = plot(NaN,NaN,'-','color',linecolors(16,:),'linewidth',2);
                    h2(10) = plot(NaN,NaN,'-','color',linecolors(18,:),'linewidth',2);
                    h2(11) = plot(NaN,NaN,'-','color',linecolors(20,:),'linewidth',2);
                    h2(12) = plot(NaN,NaN,'-','color',linecolors(22,:),'linewidth',2);
                    h2(13) = plot(NaN,NaN,'-','color',linecolors(24,:),'linewidth',2);
                    h2(14) = plot(NaN,NaN,'-','color',linecolors(26,:),'linewidth',2);
                    h2(15) = plot(NaN,NaN,'-','color',linecolors(28,:),'linewidth',2);
                    h2(16) = plot(NaN,NaN,'*','color',linecolors(28,:));  
                    legend(h2,'.5 km','1 km','2 km','3 km','4 km','5 km','6 km','7 km','8 km',...
                        '9 km','10 km','11 km','12 km','13 km','14 km','full measured fault length');
                end
            end

                xlabel('distance along transect A-A'' (km)')
                ylabel('RMS deviation of windows (m)')
                set(gca,'fontsize',16,'fontweight','bold') %,'yscale','log','xscale','log')
                xlim([min(DIST) max(DIST)])
                ylim([0 0.7])
                ax = gca;
                ax.XTick = [0:2:14 14.92];
                ax.XTickLabel = {0 2 4 6 8 10 12 14 14.9};
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

    elseif strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
        limit1 = 10000000;
        limit2 = [0 1]; % [.5 10]; 
        if limit2(1) == 0 && limit2(2) == 0
        else
            [W,W,W,slope_limited,x_RMS_limited,y_RMS_fit_limited] = fit_lin2log(log10(win_size),log10(RMS_win),limit2);
            H_RMSslope_limited = slope_limited;
            D_RMSslope_limited = 2 - slope_limited

        end
            [W,W,W,slope_allpoints,x_RMS_allpoints,y_RMS_fit_allpoints] = fit_lin2log(log10(win_size),log10(RMS_win),limit1);
            H_RMSslope_allpoints = slope_allpoints
            D_RMSslope_allpoints = 2- slope_allpoints
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
            D_RMSslope_limited2 = 2 - slope_limited
    elseif strcmpi(EQ,'Palu_long_disp.txt') == 1
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
    elseif strcmpi(EQ,'Baluchistan_stack_out_nutm41.txt') == 1        
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
    elseif isempty(regexp(EQ,'EMC')) ~= 1
        limit1 = 10000000;
        limit2 = [0 14];
            [W,W,W,slope_allpoints,x_RMS_allpoints,y_RMS_fit_allpoints] = fit_lin2log(log10(win_size),log10(RMS_win),limit1);
            [W,W,W,slope_limited,x_RMS_limited,y_RMS_fit_limited] = fit_lin2log(log10(win_size),log10(RMS_win),limit2);

            H_RMSslope_allpoints = slope_allpoints
            H_RMSslope_limited = slope_limited
            D_RMSslope_allpoints = 2- slope_allpoints
            D_RMSslope_limited = 2 - slope_limited
    end

%% SAVE DATA
    if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        save('Balochistan_plot_data.mat','win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','x_RMS_limited2','y_RMS_fit_limited2',...
            'H_RMSslope_allpoints','H_RMSslope_limited','H_RMSslope_limited2')
    elseif strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
        SV = ['Kumamoto_plot_data_FAULT',num2str(KUM_FAULT),'.mat'];
        save(SV,'win_size','RMS_win','DIST','sigma_RMS_dt',...
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
    elseif strcmpi(EQ,'Palu_long_disp.txt') == 1
        save('Palu_plot_data.mat','win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','H_RMSslope_allpoints',...
            'H_RMSslope_limited')
    elseif strcmpi(EQ,'Baluchistan_stack_out_nutm41.txt') == 1
        save('Balu_Zinke_plot_data.mat','win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','x_RMS_limited2','y_RMS_fit_limited2',...
            'H_RMSslope_allpoints','H_RMSslope_limited','H_RMSslope_limited2')    
    elseif isempty(regexp(EQ,'EMC')) ~= 1
        svnm = ['EMC_fault',EQ(4),'_plot_data.mat'];
        save(svnm,'win_size','RMS_win','distalong_transect_m','sigma_RMS_dt',...
            'RMS_win_sigma','x_RMS_limited','y_RMS_fit_limited','H_RMSslope_allpoints',...
            'H_RMSslope_limited')
    end

%% FIGURE 1 - PLOT STRAIN AND STRAIN GRADIENT ALONG RUPTURE
            
    % define sliding windows and plot
    strain_colors = [.9 .9 .9; .65 .65 .65; 0 0 0];
    grid_color = [0.85 0.85 0.85];
    figure(1)
    clf
        subplot(3,1,1)
            plot([0 max(distalong_transect_km)],[0 0],'k-')
            hold on
            plot([0 max(distalong_transect_km)],[1 1],'color',grid_color)
            hold on
            plot([0 max(distalong_transect_km)],[2 2],'color',grid_color)
            hold on
            plot([0 max(distalong_transect_km)],[3 3],'color',grid_color)
            hold on
            plot(distalong_transect_km,cumdisp_m,'b-','linewidth',2)
%             hold on
%             scatter(Field.F_distalongfault/1000,Field.f_Rslip/100,'r','filled')
                xlim([0 max(distalong_transect_m)/1000])
                if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1; 
                    hold on
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    ylim([-.5 4])
                else
                    hold on
                    plot([2.25 2.25],[-10 108],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    ylim([0 3.5])
                end
                ylabel('lateral fault offset (m)')
                set(gca,'fontsize',16,'fontweight','bold','layer','top','XTickLabel',[])
        subplot(3,1,2)
            for pp = 1: length(define_win)
                plot(Slide(pp).x_strain/1000,Slide(pp).strain,'color',strain_colors(pp,:),'linewidth',2)
                hold on
            end
%             plot([0 max(distalong_transect_m)/1000],[3.2e-04 3.2e-04],'r--')
                xlim([0 max(distalong_transect_m)/1000])
                if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1; 
                    plot([0 max(Slide(1).x_strain)],[(mean([Slide(3).strain])+std([Slide(3).strain])) (mean([Slide(3).strain])+std([Slide(3).strain]))],'r--')
                    hold on
                    plot([0 max(Slide(1).x_strain)],[(mean([Slide(3).strain])-std([Slide(3).strain])) (mean([Slide(3).strain])-std([Slide(3).strain]))],'r--')
                    hold on
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    ylim([-0.02 0.02]);                
                else
                    hold on
                    plot([2.25 2.25],[-10 108],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    ylim([-.02 .02])
                end
%                 legend('3','5','10','1 standard deviation, 10 point window')
                ylabel('strain')
                set(gca,'fontsize',16,'fontweight','bold','layer','top','XTickLabel',[])
        subplot(3,1,3)
            for qq = 1: length(define_win)
                plot(Slide(qq).x_sgrad/1000,Slide(qq).strain_grad,'color',strain_colors(qq,:),'linewidth',2)
                hold on
            end
                xlim([0 max(distalong_transect_m)/1000])
                if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1; 
                    plot([0 max(Slide(1).x_sgrad)/1000],[(mean([Slide(3).strain_grad])+std([Slide(3).strain_grad])) (mean([Slide(3).strain_grad])+std([Slide(3).strain_grad]))],'r--')
                    hold on
                    plot([0 max(Slide(1).x_sgrad)/1000],[(mean([Slide(3).strain_grad])-std([Slide(3).strain_grad])) (mean([Slide(3).strain_grad])-std([Slide(3).strain_grad]))],'r--')
                    hold on
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    ylim([-5e-4 5e-4])
                else
                    hold on
                    plot([2.25 2.25],[-10 108],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    ylim([-5e-4 5e-4])
                end
                legend('3 point average','5 point average','10 point average','1 standard deviation, 10 point average')
                xlabel('distance along transect A-A'' (km)')
                ylabel('strain gradient (m^{-1})')
                set(gca,'fontsize',16,'fontweight','bold','layer','top')
             
%% FIGURE 2 - PSD THOMSON
    figure (2)
    clf
    subplot(2,3,1:3)
        % plot grid lines
            plot([0 max(distalong_transect_km)],[-1 -1],'color',grid_color)
            hold on
            plot([0 max(distalong_transect_km)],[0 0],'k-')
            hold on
            plot([0 max(distalong_transect_km)],[1 1],'color',grid_color)
            hold on
            plot([0 max(distalong_transect_km)],[2 2],'color',grid_color)
            hold on
            plot([0 max(distalong_transect_km)],[3 3],'color',grid_color)
            hold on
        % plot transect bends
        if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
            plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
            hold on
            plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
            hold on
        end
        % plot data
            scatter(distalong_transect_km,cumdisp_m,'b','filled')
            hold on
            if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
                errorbar(distalong_transect_km,cumdisp_m,err,'b')
                hold on
            end
            scatter(f_dist_km,f_Rslip/100,'r','filled')        
                ylabel('lateral fault offset (m)')
                set(gca,'fontsize',16,'fontweight','bold')
                xlim([0 max(distalong_transect_km)])
                xlabel('distance along transect A-A'' (km)')
                ax = gca;
                ax.XTick = [0:2:14 14.92];
                ax.XTickLabel = {0 2 4 6 8 10 12 14 14.9};
                ylim([-0.5 3.5])
                xlim([0 max(distalong_transect_km)])
    subplot(2,3,4)
        xt = 10^4.2;
        yt = 10^-1;
        if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
            loglog(1./f,p,'k');
            hold on
            xlabel('wavelength (m)')
            set(gca,'fontsize',16,'fontweight','bold')
        elseif strcmpi(PSD_xaxis,'frequency') == 1;
            loglog(f,p,'k');
            hold on
            xlabel('wavenumber (m^{-1})')
            set(gca,'xdir','reverse','fontsize',16,'fontweight','bold')
        end
            loglog(10.^x_fft,10.^yCalc_fft,'r-','linewidth',2)

            text_1 = ['D = ',num2str(D_fft)];
            text_2 = ['H = ',num2str(H_fft)];
            text_3 = ['R^2 = ',num2str(RSq_fft)];
            text(xt,yt,{text_1,text_2,text_3},'fontsize',16,'fontweight','bold')

            ylabel('power spectral density (m^3)')
            title('Peridiogram method') 
            set(gca,'fontsize',16,'fontweight','bold')
%             legend('PSD','best fit')
            xlim([60 20000])
            ylim([10^-2 10^4])
            grid on
            
    subplot(2,3,5)
        xt = 10^4.2;
        yt = 10^0;
        if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
            loglog(1./freq,(pxx),'k')
            hold on
            xlabel('wavelength (m)')
            set(gca,'fontsize',16,'fontweight','bold')
        elseif strcmpi(PSD_xaxis,'frequency') == 1;
            loglog(freq,(pxx),'k')
            hold on
            xlabel('wavenumber (m^{-1})')
            set(gca,'xdir','reverse','fontsize',16,'fontweight','bold')
        end    
            loglog(10.^x_Thom,10.^yCalc_Thom,'r-','linewidth',2)
    % 
            text_1 = ['D = ',num2str(D_Thom)];
            text_2 = ['H = ',num2str(H_Thom)];
            text_3 = ['R^2 = ',num2str(RSq_Thom)];
            text(xt,yt,{text_1,text_2,text_3},'fontsize',16,'fontweight','bold')

%             ylabel('power spectral density (m^3)')
    %         title({'power spectral density - Thomson Multitaper';text_1;text_2}) 
            title({'Thomson multitaper method'})
            set(gca,'fontsize',16,'fontweight','bold')
%             legend('PSD','best fit')
            xlim([60 20000])
            ylim([10^-2 10^4])
            grid on
     
     subplot(2,3,6)       
        xt = 10^4.2;
        yt = 10^0;
        if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
            loglog(1./Fw,Pyy,'k');
            hold on
            xlabel('wavelength (m)')
            set(gca,'fontsize',16,'fontweight','bold')
        elseif strcmpi(PSD_xaxis,'frequency') == 1;
            loglog(Fw,Pyy,'k');
            hold on
            xlabel('wavenumber (m^{-1})')
            set(gca,'xdir','reverse','fontsize',16,'fontweight','bold')        
        end
            loglog(10.^x_welch,10.^yCalc_welch,'r-','linewidth',2)

            text_1 = ['D = ',num2str(D_welch)];
            text_2 = ['H = ',num2str(H_welch)];
            text_3 = ['R^2 = ',num2str(RSq_welch)];
            text(xt,yt,{text_1,text_2,text_3},'fontsize',16,'fontweight','bold')

%             ylabel('power spectral density (m^3)')
            title('Welch method') 
            xlim([60 20000])
            ylim([10^-2 10^4])
%             legend('PSD','best fit')
            grid on
            
%% FIGURE 3 - PSD ALL
    figure (3)
    clf
        if strcmpi(PSD_xaxis,'wavelength') == 1; % "wavelength" or "frequency"
            loglog(1./f,p,'r-','linewidth',2);
            hold on
            loglog(1./freq,(pxx),'r-','linewidth',2)
            hold on
            loglog(1./Fw,Pyy,'r-','linewidth',2);
            hold on
            xlabel('wavelength (m)')
            set(gca,'fontsize',12,'fontweight','bold')
        elseif strcmpi(PSD_xaxis,'frequency') == 1;
            loglog(f,p,'r-','linewidth',2);
            hold on
            loglog(freq,(pxx),'r-','linewidth',2)
            hold on
            loglog(Fw,Pyy,'r-','linewidth',2);
            hold on
            xlabel('wavenumber (m^{-1})')
            set(gca,'xdir','reverse','fontsize',12,'fontweight','bold')
        end
            loglog(10.^x_fft,10.^yCalc_fft,'r-','linewidth',2)
            hold on
            loglog(10.^x_Thom,10.^yCalc_Thom,'r-','linewidth',2)
            hold on
            loglog(10.^x_welch,10.^yCalc_welch,'r-','linewidth',2)
            
            ylabel('power spectral density (m^3)')
            grid on
            xlim([60 20000])
            
%% FIGURE 4 - RMS,D, H v length - log
    figure(4)
    clf
    plot(win_size,H_fft_win,'b-','linewidth',2)
    hold on
    plot(max(distalong_transect_m)/1000,H_fft,'b*','MarkerSize',12)
    hold on
    plot(win_size,D_fft_win,'r-','linewidth',2)
    hold on
    plot(max(distalong_transect_m)/1000,D_fft,'r*','MarkerSize',12)
    hold on
    plot(win_size,RMS_win,'g-','linewidth',2)
    hold on
    plot(max(distalong_transect_m)/1000,sigma_RMS_dt,'g*','MarkerSize',12)
    hold on
    errorbar(win_size,RMS_win,RMS_win_sigma,'color','g')
    
    if strcmpi(EQ,'Balochistan_Vallage_FarField.txt') == 1
        xlim([4 200])
        ylim([0.02 3.1])
%         set(gca, 'XTick',[0:15]); % 0 1 2 3 4 5 6 7 8 9 10 12 15])
%         set(gca, 'YTick',[0:.1:1])
    elseif strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
        hold on
        plot(40/1000,0.23,'go')
        xlim([40/1000 (max(distalong_transect_m)/1000+1)])
        ylim([-1 2.5])
%         set(gca, 'XTick',[0:15]); % 0 1 2 3 4 5 6 7 8 9 10 12 15])
%         set(gca, 'YTick',[0:.1:1])
    elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
        xlim([.4 (max(distalong_transect_m)/1000+1)])
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
    figure(6)
    clf
    plot(win_size,H_fft_win,'b-','linewidth',2)
    hold on
    plot(max(distalong_transect_m)/1000,H_fft,'b*','MarkerSize',12)
    hold on
    plot(win_size,D_fft_win,'r-','linewidth',2)
    hold on
    plot(max(distalong_transect_m)/1000,D_fft,'r*','MarkerSize',12)
    hold on
    plot(win_size,RMS_win,'g-','linewidth',2)
    hold on
    plot(max(distalong_transect_m)/1000,sigma_RMS_dt,'g*','MarkerSize',12)
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
        ax = gca;
        ax.XTick = [0:2:14 14.92];
        ax.XTickLabel = {0 2 4 6 8 10 12 14 14.9};
        grid on
        
%% FIGURE 7 - HURST FROM RMS ROUGHNESS
    figure(7)
    clf
    plot(win_size,RMS_win,'g-','linewidth',2)
    hold on
    plot(max(distalong_transect_m)/1000,sigma_RMS_dt,'g*','MarkerSize',12)
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
    elseif strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
        hold on
        plot(40/1000,0.23,'go')
        xlim([40/1000 (max(distalong_transect_m)/1000+1)])
        ylim([-1.5 10^-.1])
%         set(gca, 'XTick',[0:15]); % 0 1 2 3 4 5 6 7 8 9 10 12 15])
%         set(gca, 'YTick',[0:.1:1])
    elseif strcmpi(EQ,'Landers_disp_data_no_header.csv') == 1
        xlim([.4 (max(distalong_transect_m)/1000+1)])
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
    figure(8)
    clf
    % plot grid lines
        plot([0 max(distalong_transect_km)],[-1 -1],'color',grid_color)
        hold on
        plot([0 max(distalong_transect_km)],[0 0],'k-')
        hold on
        plot([0 max(distalong_transect_km)],[1 1],'color',grid_color)
        hold on
        plot([0 max(distalong_transect_km)],[2 2],'color',grid_color)
        hold on
        plot([0 max(distalong_transect_km)],[3 3],'color',grid_color)
        hold on
    % plot transect bends
        if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
            plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
            hold on
            plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
            hold on
        end
    % plot data
        plot(DIST,DISP,'b-','linewidth',2)
        hold on
        plot(DIST,profile,'g-','linewidth',2)
        hold on
        plot(DIST,profile_taper,'r-','linewidth',2)
        hold on
%             tit_text = ['Slip distributions for ',EQ_name];
%             title(tit_text)
            xlim([0 max(DIST)])
            ylim([-2 3.5])
            xlabel('distance along transect A-A'' (km)')
            ylabel('lateral fault offset (m)')
            set(gca,'fontsize',16,'fontweight','bold')
            ax = gca;
            ax.XTick = [0:2:14 14.92];
            ax.XTickLabel = {0 2 4 6 8 10 12 14 14.9};
            
            h8 = zeros(3, 1);
            h8(1) = plot(NaN,NaN,'b-','linewidth',2);
            h8(2) = plot(NaN,NaN,'g-','linewidth',2);
            h8(3) = plot(NaN,NaN,'r-','linewidth',2);
            legend(h8,'raw slip','detrended slip','detrended and tapered slip')
        
%% FIGURE 9 - plot cumulative displacement with RMS and mean slip

if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
    if KUM_FAULT == 1 || KUM_FAULT == 2 || KUM_FAULT == 3 || ...
                        KUM_FAULT == 4 || KUM_FAULT == 5 || KUM_FAULT == 6;
    else
        fig = figure(9);
        clf
        set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
        yyaxis left
                % plot data
                    plot(distalong_transect_km,cum_Disp,'-','color',[.8 .8 .8],'linewidth',1)
                    hold on
                    plot(win_center_1km,mean_slip_1km,'r--','linewidth',2)
                    hold on
                    plot(win_center_3km,mean_slip_3km,'c--','linewidth',2)
                    hold on
                    plot(win_center_5km,mean_slip_5km,'b--','linewidth',2)
               % plot grid lines
                    plot([0 max(distalong_transect_km)],[0 0],'k-')
                    hold on
                    plot([0 max(distalong_transect_km)],[1 1],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[2 2],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[3 3],'color',grid_color)
                    hold on
               % plot transect bends
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                        ylim([-.5 3.5])
                        ylabel('lateral fault offset (m)')
        yyaxis right
            plot(win_center_1km,RMS_1km,'r-','linewidth',2)
            hold on
            plot(win_center_3km,RMS_3km,'c-','linewidth',2)
            hold on
            plot(win_center_5km,RMS_5km,'b-','linewidth',2)
            hold on 
            plot(NaN,NaN,'r*')
            hold on 
            plot(NaN,NaN,'c*')
            hold on 
            plot(NaN,NaN,'b*')

            ylabel('RMS deviation (m)')
            text_1km = ['R = ',num2str(R_RMSdisp(2)),' P = ',num2str(P_RMSdisp(2))];
            text_3km = ['R = ',num2str(R_RMSdisp(6)),' P = ',num2str(P_RMSdisp(6))];
            text_5km = ['R = ',num2str(R_RMSdisp(10)),' P = ',num2str(P_RMSdisp(10))];
            xlim([0 max(distalong_transect_km)])
            set(gca,'fontsize',20,'fontweight','bold')
            xlabel('distance along transect A-A'' (km)')
            ax = gca;
            ax.XTick = [0:2:14 14.92];
            ax.XTickLabel = {0 2 4 6 8 10 12 14 14.9};

            legend('measured slip distribution','mean slip 1 km window','mean slip 3 km window',...
                'mean slip 5 km window','RMS deviation 1 km window','RMS deviation 3 km window',...
                'RMS deviation 5 km window',text_1km,text_3km,text_5km)
    end
end

%% FIGURE 10 - plot cumulative displacement with RMS and mean slip
if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
    if KUM_FAULT == 1 || KUM_FAULT == 2 || KUM_FAULT == 3 || ...
                        KUM_FAULT == 4 || KUM_FAULT == 5 || KUM_FAULT == 6;
    else

        color1 = [254 196 79]/256; %'r'; 
        color2 = [0 136 55]/256; %'b'; 
        color3 = [123 50 148]/256; %'c';
        figure(10);
        clf
        fig1 = subplot(2,1,1);
            yyaxis left
                % plot grid lines
                    plot([0 max(distalong_transect_km)],[0 0],'k-')
                    hold on
                    plot([0 max(distalong_transect_km)],[1 1],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[2 2],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[3 3],'color',grid_color)
                    hold on
                % plot transect bends
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                % plot data
                    plot(distalong_transect_km,cum_Disp,'-','color',[.8 .8 .8],'linewidth',1)
                    hold on
                    plot(win_center_1km,mean_slip_1km,'-','color',color1,'linewidth',2)
                    hold on
                    plot(win_center_3km,mean_slip_3km,'--','color',color1,'linewidth',2)
                    hold on
                    plot(win_center_5km,mean_slip_5km,':','color',color1,'linewidth',2)
                        ax11 = gca;
                        ax11.YColor = [0 0 0];
                        ylim([-.5 4])
                        ylabel('lateral fault offset (m)')
            yyaxis right
                plot(win_center_1km,mean_str_1km,'-','color',color2,'linewidth',2)
                hold on
                plot(win_center_3km,mean_str_3km,'--','color',color2,'linewidth',2)
                hold on
                plot(win_center_5km,mean_str_5km,':','color',color2,'linewidth',2)
                hold on
                    ax1 = gca;
                    ax1.XTick = ([0:2:15 15.48]);
                    ax1.XTickLabel = [];
                    ax1.YTick = ([0:0.004:0.016]);
                    ax1.YColor = [0 0 0];
                    ylabel('strain magnitude')
                    xlim([0 max(distalong_transect_km)])
                    ylim([0 0.016])
                    set(gca,'fontsize',16,'fontweight','bold')

        fig2 = subplot(2,1,2);
            yyaxis left
                % plot grid lines
                    plot([0 max(distalong_transect_km)],[0 0],'k-')
                    hold on
                    plot([0 max(distalong_transect_km)],[1 1],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[2 2],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[3 3],'color',grid_color)
                    hold on
                % plot transect bends
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                % plot data
                    plot(distalong_transect_km,cum_Disp,'-','color',[.8 .8 .8],'linewidth',1)
                    hold on
                    plot(win_center_1km,mean_slip_1km,'-','color',color1,'linewidth',2)
                    hold on
                    plot(win_center_3km,mean_slip_3km,'--','color',color1,'linewidth',2)
                    hold on
                    plot(win_center_5km,mean_slip_5km,':','color',color1,'linewidth',2)
                        ax21 = gca;
                        ax21.YColor = [0 0 0];
                        ylim([-.5 4])
                        ylabel('lateral fault offset (m)')
            yyaxis right
                    plot(win_center_1km,RMS_1km,'-','color',color3,'linewidth',2)
                    hold on
                    plot(win_center_3km,RMS_3km,'--','color',color3,'linewidth',2)
                    hold on
                    plot(win_center_5km,RMS_5km,':','color',color3,'linewidth',2)
                        ax2 = gca;
                        ax2.XTick = [0:2:14 14.92];
                        ax2.XTickLabel = {0 2 4 6 8 10 12 14 14.9};
                        ax2.YColor = [0 0 0];
                        ylabel('RMS deviation (m)')
                        xlim([0 max(distalong_transect_km)])
                        set(gca,'fontsize',16,'fontweight','bold')
                        xlabel('distance along transect A-A'' (km)')

                        h2 = zeros(3, 1);
                        h2(1) = plot(NaN,NaN,'-','color',[.8 .8 .8],'linewidth',1);
                        h2(2) = plot(NaN,NaN,'-','color',color1,'linewidth',2);
                        h2(3) = plot(NaN,NaN,'-','color',color2,'linewidth',2);
                        h2(4) = plot(NaN,NaN,'-','color',color3,'linewidth',2);
                        legend(h2, 'measured fault offset','average fault offset','average strain','average \sigma_{RMS}');

                display(['Offset - Strain 1km: R = ',num2str(R_StrDisp(2)),' P = ',num2str(P_StrDisp(2))]);
                display(['Offset - Strain 3km: R = ',num2str(R_StrDisp(6)),' P = ',num2str(P_StrDisp(6))]);
                display(['Offset - Strain 5km: R = ',num2str(R_StrDisp(10)),' P = ',num2str(P_StrDisp(10))]);
                display(['Offset - RMS 1km: R = ',num2str(R_RMSdisp(2)),' P = ',num2str(P_RMSdisp(2))]);
                display(['Offset - RMS 3km: R = ',num2str(R_RMSdisp(6)),' P = ',num2str(P_RMSdisp(6))]);
                display(['Offset - RMS 5km: R = ',num2str(R_RMSdisp(10)),' P = ',num2str(P_RMSdisp(10))]);
                display(['RMS - Strain 1km: R = ',num2str(R_StrRMS(2)),' P = ',num2str(P_StrRMS(2))]);
                display(['RMS - Strain 3km: R = ',num2str(R_StrRMS(6)),' P = ',num2str(P_StrRMS(6))]);
                display(['RMS - Strain 5km: R = ',num2str(R_StrRMS(10)),' P = ',num2str(P_StrRMS(10))]);
    end
end

%% FIGURE 11 - PLOT MEASUREMENT LOCATIONS
if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
    
    if KUM_FAULT == 1 || KUM_FAULT == 2 || KUM_FAULT == 3 || ...
                        KUM_FAULT == 4 || KUM_FAULT == 5 || KUM_FAULT == 6;
    else     
        dot_sz = 20
        figure(11)
        clf
            % plot transect A - A'
                plot(prof_x,prof_y,'--','color',[0.6 0.6 0.6],'linewidth',2)
                hold on
            % plot field measurement locations
                scatter(Field.f_easting,Field.f_northing,10,'markerfacecolor','r','markeredgecolor','k')
            % plot ICP measurement locations
            for ww = 1:max(size(Faults))
                scatter(Scarps_x(Faults(ww).ind),Scarps_y(Faults(ww).ind),dot_sz,'markerfacecolor',Faults(ww).color,'markeredgecolor','none')
                hold on
            end

                set(gca,'fontsize',18,'fontweight','bold')
                xlabel('easting (m)')
                ylabel('northing (m)')
                ax = gca;
                ax.XTick = [-20000 -18000 -16000 -14000 -12000 -10000 -8000];
                ax.XTickLabel = {-20000 -18000 -16000 -14000 -12000 -10000 -8000};
                ax.XTickLabelRotation = 45;
                ax.YTick = [-30000 -28000 -26000 -24000 -22000 -20000];
                ax.YTickLabel = {-30000 -28000 -26000 -24000 -22000 -20000};
                ax.YTickLabelRotation = 45;
                ax.YMinorTick = 'on';
                ax.XMinorTick = 'on';
                axis equal
                xlim([min(prof_x)-1000 max(prof_x)+1000])
                ylim([min(prof_y)-1000 max(prof_y)+1000])
    end
end

%% FIGURE 12 - INDIVIDUAL FAULTS OFFSET
if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
    if KUM_FAULT == 1 || KUM_FAULT == 2 || KUM_FAULT == 3 || ...
                        KUM_FAULT == 4 || KUM_FAULT == 5 || KUM_FAULT == 6;
    else
        pt_sz = 10;
        ln_wd = 2;
        figure (12)
        clf
                plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                hold on
                plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                hold on
                plot([0 max(distalong_transect_km)],[0 0],'k-')
                hold on
                plot([0 max(distalong_transect_km)],[-1 -1],'color',[.8 .8 .8])
                hold on
                plot([0 max(distalong_transect_km)],[1 1],'color',[.8 .8 .8])
                hold on
                plot([0 max(distalong_transect_km)],[2 2],'color',[.8 .8 .8])
                hold on
                plot([0 max(distalong_transect_km)],[3 3],'color',[.8 .8 .8])
                hold on
                plot(distalong_transect_km,cumdisp_m,'color',[.6 .6 .6])
                hold on
                scatter(f_dist_km,f_Rslip/100,'r','filled')
                hold on
                for xx = 1:max(size(Faults))
                    plot(Scarps_distalongtransect(Faults(xx).ind),Scarps_lat(Faults(xx).ind),'color',Faults(xx).color,'linewidth',ln_wd)
                    hold on
                end

                    xlim([0 max(distalong_transect_km)])
                    ylim([-1.5 3.5])
                    ax = gca;
                    ax.XTick = [0:2:14 14.92];
                    ax.XTickLabel = {0 2 4 6 8 10 12 14 14.9};
                    set(gca,'fontsize',18,'fontweight','bold')
                    xlabel('distance along transect A-A'' (km)')
                    ylabel('lateral fault offset (m)')
    end
end
        
%% FIGURE 13 - INDIVIDUAL FAULTS STRAIN
if strcmpi(EQ,'Kumamoto_40m_3seg_AFTERBUG.mat') == 1
    if KUM_FAULT == 1 || KUM_FAULT == 2 || KUM_FAULT == 3 || ...
                        KUM_FAULT == 4 || KUM_FAULT == 5 || KUM_FAULT == 6;
    else
               % define sliding windows and plot
        grid_color = [0.85 0.85 0.85];
        WF_color = [0.8 0.8 0.8];
        ln_wd = 2;
        figure(13)
        clf
            subplot(4,1,1:2)
                % plot grid lines
                    plot([0 max(distalong_transect_km)],[-1 -1],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[0 0],'k-')
                    hold on
                    plot([0 max(distalong_transect_km)],[1 1],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[2 2],'color',grid_color)
                    hold on
                    plot([0 max(distalong_transect_km)],[3 3],'color',grid_color)
                    hold on
                % plot fault segment breaks
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                % plot field data
                    scatter(f_dist_km,f_Rslip/100,'r','filled')
                    hold on
                % plot whole fault offset
                    plot(distalong_transect_km,cumdisp_m,'-','color',WF_color,'linewidth',1)
                    hold on
                    for xx = 1:max(size(Faults))
                        plot(Scarps_distalongtransect(Faults(xx).ind),Scarps_lat(Faults(xx).ind),'color',Faults(xx).color,'linewidth',ln_wd)
                        hold on
                    end

                        ylim([-1.5 3.5])
                        xlim([0 max(distalong_transect_m)/1000])
                        ylabel('lateral fault offset (m)')
                        set(gca,'fontsize',16,'fontweight','bold','layer','top','XTickLabel',[])
            subplot(4,1,3)
                % plot zero line
                    plot([0 max(distalong_transect_km)],[0 0],'color','k')
                    hold on
                % plot fault segment breaks
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                % plot standard deviation of strain for whole fault
                    plot([0 max(Slide(1).x_strain)],[(mean([Slide(strain_plot).strain])+std([Slide(strain_plot).strain])) (mean([Slide(strain_plot).strain])+std([Slide(strain_plot).strain]))],'r--')
                    hold on
                    plot([0 max(Slide(1).x_strain)],[(mean([Slide(strain_plot).strain])-std([Slide(strain_plot).strain])) (mean([Slide(strain_plot).strain])-std([Slide(strain_plot).strain]))],'r--')
                    hold on
                % plot strains
                    plot(Slide(strain_plot).x_strain/1000,Slide(strain_plot).strain,'color',WF_color,'linewidth',1)
                    hold on
                    for xx2 = 1:max(size(Faults))
                        plot(Faults(xx2).Slide(strain_plot).x_strain,Faults(xx2).Slide(strain_plot).strain,'color',Faults(xx2).color,'linewidth',ln_wd)
                        hold on
                    end
                        if strain_plot == 3
                            ylim([-7e-3 7e-3]); 
                        elseif strain_plot == 2
                            ylim([-1.2e-2 1.2e-2]); 
                        end
                        xlim([0 max(distalong_transect_m)/1000])
                        ylabel('strain')
                        set(gca,'fontsize',16,'fontweight','bold','layer','top','XTickLabel',[])
            subplot(4,1,4)
                % plot zero line
                    plot([0 max(distalong_transect_km)],[0 0],'color','k')
                    hold on
                % plot fault segment breaks
                    plot([trans_bend1 trans_bend1],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                    plot([trans_bend2 trans_bend2],[-10 10],'linestyle','--','color',[0.6 0.6 0.6],'linewidth',2)
                    hold on
                % plot standard deviation of strain gradient for whole fault
                    plot([0 max(Slide(1).x_sgrad)/1000],[(mean([Slide(strain_plot).strain_grad])+std([Slide(strain_plot).strain_grad])) (mean([Slide(strain_plot).strain_grad])+std([Slide(strain_plot).strain_grad]))],'r--')
                    hold on
                    plot([0 max(Slide(1).x_sgrad)/1000],[(mean([Slide(strain_plot).strain_grad])-std([Slide(strain_plot).strain_grad])) (mean([Slide(strain_plot).strain_grad])-std([Slide(strain_plot).strain_grad]))],'r--')
                    hold on
                % plot strain gradientts
                    plot(Slide(strain_plot).x_sgrad/1000,Slide(strain_plot).strain_grad,'color',WF_color,'linewidth',1)
                    hold on
                    for xx2 = 1:max(size(Faults))
                        plot(Faults(xx2).Slide(strain_plot).x_sgrad,Faults(xx2).Slide(strain_plot).strain_grad,'color',Faults(xx2).color,'linewidth',ln_wd)
                        hold on
                    end
                        if strain_plot == 3
                            ylim([-2e-4 2e-4]); 
                        elseif strain_plot == 2
                            ylim([-3.5e-4 3.5e-4]); 
                        end
                        xlim([0 max(distalong_transect_m)/1000])
    %                 legend('3 point average','5 point average','10 point average','1 standard deviation, 10 point average')
                    ax = gca;
                    ax.XTick = [0:2:14 14.92];
                    ax.XTickLabel = {0 2 4 6 8 10 12 14 14.9};
                    xlabel('distance along transect A-A'' (km)')
                    ylabel('strain gradient (m^{-1})')
                    set(gca,'fontsize',16,'fontweight','bold','layer','top')    
    end       
    end

