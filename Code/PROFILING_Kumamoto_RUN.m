% Lia Lajoie
% Colorado School of Mines
% 17 Aug 2016

%% ************************************************************************               
% Inputs are .csv files
% with easting (m) and northing (m) values for fault traces (each fault
% section has own file with all x and y values) and structure storing x and
% y limits of ICP data to plot as well as the reference section numbers
% from Fletcher fault traces.
% *************************************************************************
% link to other functions
addpath('C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\ElMayorCucapah');
addpath('C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Scripts');
% shut down all figuresI
fh=findall(0,'type','figure');
for i=1:length(fh)
     clo(fh(i));
     clf
end

clear all

%% PROFILE TO START AT 
% profiles advance along transect from SW to NE
% use this if the first profiles in dataset do not have meaningful data and
% you want to skip over them

% load Kumamoto_40width
PROF_START = 1; 

%% USER DEFINES FAULT TO RUN AND METHODS
FAULT = 11; % which fault segment to analyze 10 = all faults, 11 = generalized fault with 3 segs
PROF_ANALYZE = 'off'; % turn profile analyzer 'on' or 'off'
    prof_len = 1500; % half length of profiles (in m)
    prof_width = 20; % half width of profiles (in m)
PROF_MODE = 'even'; % data = center on data; even = evenly spaced
    prof_space = prof_width*2; % profile spacing, meters - only used for PROF_MODE = 'even'
PROFDATA_MODE = 'prof_determined'; % global = data rotated to GLOBAL fault azimuth; data_determined = data rotated to closest fault segment to DATA POINT; prof_determined = data rotated to closest fault segment to PROFILE
BIN_DIST = 'off'; % generate plots with discrete distance bins?
    bin_size = 100; % width of distance-from-fault bins
    num_bins = 10; 
    smooth_win = 20; % smoothing window for along fault profile curves
DATA_MODE = 'local'; % global = data parallel to whole fault; local = data parallel to local fault
    global_fault_az = 45; % global azimuth of main fault - only used for DATA_MODE = 'global'
   
% open and extract fault name data for FAULT selected above
FaultData_Kumamoto_store;
NAME = strcat(FaultData(FAULT).FaultName,' Fault; ',FaultData(FAULT).SegmentName,' Segment');

%% DEFINE VARIABLES 
gridsize = 8; % (m) size of pixel tracking grid 
n = 1000; % max allowable picks per fault segment

%% ************************************************************************
%                    IMPORT LAJOIE COSI-CORR DATA
%**************************************************************************
% Load COSI data for specified fault
    folderpath = 'C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Cosi-Corr\FaultSegments\';
    filename = ['COSIdata_fault',num2str(FAULT),'.mat'];
    loadfile = [folderpath,filename];
    load(loadfile);
            % NS_grid/ EW_grid: NS/ EW displacements computed by COSI-Corr
            % SNR_grid: signal to noise ratio 
            % x_grid, y_grid
            % x_vect, y_vect

%% ************************************************************************
%                          IMPORT FIELD DATA
%**************************************************************************
% SPECIFIC DATA POINTS TO FAULT (JUST X,Y)
field_folderpath = 'C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Data\FieldData\FAULT_segments';
filestart = 'FieldData_Fault';
fileend = '.csv';
    file_name = [field_folderpath,'\',filestart,num2str(FAULT),fileend]; % get file path name
    field_fault_data = importdata(file_name);
    x_field_fault = field_fault_data(:,1);
    y_field_fault = field_fault_data(:,2);

% FIELD DATA - ENTIRE DATASET, NOT FAULT SPECIFIC
    ParseField_completedataset;
    % Field_Rslip; Field_Sup; Field_dipslip; Field_northing; Field_easting
   
% CROSS REFERENCE FIELD DATA
% Cross reference field data. I  plotted all field measurements
% in Matlab and filtered datapoints to only those located on each specific 
% fault trace. When saving the data, only x and y locations are retained, 
% so this step cross-references the x,y locations of field data on each 
% fault to the full dataset and extracts other data values.
    field_indeces = [];
    for t = 1:length(x_field_fault)
        x_find = find(Field_easting == x_field_fault(t) & Field_northing == y_field_fault(t));
        field_indeces = [field_indeces x_find'];
    end
        
% COLLECT FAULT SPECIFIC FIELD DATA
    f_easting = Field_easting(field_indeces); 
    f_northing = Field_northing(field_indeces); 
    f_Rslip = Field_Rslip(field_indeces); 
    f_Sup = Field_Sup(field_indeces); 
    f_dipslip = Field_dipslip(field_indeces); 
    
%% ************************************************************************
%                   CALCULATE FAULT SEGMENT DATA
%**************************************************************************    
% create vectors of fault endpoints (x and y values along fault trace)
    Fxx = FaultData(FAULT).Xfault;
    Fyy = FaultData(FAULT).Yfault;
    numsegs_FAULT = length(Fxx)-1; % number of fault segments

% ORDER POINTS NORTH TO SOUTH - ALL CALCULATIONS ASSUME THIS POINT ORDER
    [Fy,ind] = sort(Fyy,'descend');
    Fx = Fxx(ind);
   
% make vector of start and end points for FAULT segments
    Sx_end1 = Fx(1:end-1);
    Sy_end1 = Fy(1:end-1);
    Sx_end2 = Fx(2:end);
    Sy_end2 = Fy(2:end);
  
% calculate dx, dy, scale factor (a), and length for each fault segment
    dx = Sx_end2-Sx_end1; % dx - vector with length = number of segments
    dy = Sy_end2-Sy_end1; % dy
    a_seg = (prof_space^2./(dx.^2 + dy.^2)).^(1/2); % reduction factor to convert the vector dxdy to a length of profile_separation
    seglen = sqrt((dx).^2+(dy).^2); % length of each segment
    m = dy./dx; % slope of each segment
    norm = -(m).^-1; % normal to each segment
    
    TotalSegDist = sum(seglen); 
    tot_proflen = cumsum(seglen);
    add_proflen = [0 tot_proflen];
%     add_proflen = add_proflen(1:end-1);
    numprofs = ceil(TotalSegDist/prof_space);
    CumProfDist = cumsum(ones(numprofs,1)*prof_space)-prof_space; % profile distances from starting point
    profnum = 1:numprofs;
    
%         seglen = sqrt((Sx_end1-Sx_end2).^2+(Sy_end1-Sy_end2).^2);
%         tot_proflen = cumsum(seglen); 
%         add_proflen = [0 tot_proflen'];
    
% calculate fault azimuth for each fault segment
    local_fault_az = zeros(size(a_seg));
    for ii = 1:numsegs_FAULT
        if dx(ii) < 0 
            local_fault_az(ii) = 90 - atand(dx(ii)/dy(ii));
        elseif dx(ii) > 0 
            local_fault_az(ii) = 90 + atand(dy(ii)/dx(ii));
        end
    end
     
%% ************************************************************************
%                      FIND PROFILE LOCATIONS
%**************************************************************************
    if strcmp(PROF_MODE,'even') == 1;
        numseg = numsegs_FAULT;
        profsegs = zeros(numprofs,1); 
        a_start = zeros(numseg,1);
        slopeseg = zeros(numprofs,1);
        normseg = zeros(numprofs,1);
        startprof = ones(numseg,1);
        trackseg = zeros(numseg,1);
        countprofs = trackseg;

    % compute number of profiles and key each profile to its host
    % segment
        for jj = 1:numseg
            % Original 
            if jj == 1
                profs = CumProfDist <= tot_proflen(jj); % find profiles in segment j
            elseif jj > 1
                profs = CumProfDist > tot_proflen(jj-1) & CumProfDist <= tot_proflen(jj);
            end
            profsegs(profs) = jj; % segment each profile belongs to
            slopeseg(profs) = m(jj); % slope of each profile in segment
            normseg(profs) = norm(jj); % normal of each profile in degment
            countprofs(jj) = sum(profs);

            % identify which segments have profiles
            if sum(profs) > 0
                endprof = max(profnum(profs)); % last profile in segment
                s_prof = min(profnum(profs));
                startprof(jj) = s_prof;
                startdist = CumProfDist(s_prof) - add_proflen(jj);
                a_start(jj) = sqrt(startdist^2/(dx(jj)^2 + dy(jj)^2));
                trackseg(jj) = 1;
            else
                trackseg(jj) = 0;
            end
        end

    % Find Profile Centers
        Sx = zeros(numprofs,1);
        Sy = Sx;
        for kk = 1:numprofs
            seg = profsegs(kk);
            Sx(kk) = Fx(seg) + a_start(seg)*dx(seg) + (kk-startprof(seg))*a_seg(seg)*dx(seg); % x-coordinate of profile center
            Sy(kk) = Fy(seg) + a_start(seg)*dy(seg) + (kk-startprof(seg))*a_seg(seg)*dy(seg); % y-coordinate of profile center
        end

        figure(2)
        clf
        plot(Fx,Fy)
        hold on
        plot(Sx,Sy,'ro');
        axis equal

    %         Sx = [XX1:a_seg*dx:XX2]';
    %         Sy = [YY1:a_seg*dy:YY2]';

elseif strcmp(PROF_MODE,'data') == 1;      
    Sx = f_easting;
    Sy = f_northing;
end
    
%     % Segment end points (between profile locations)
%         Sx_end1 = Sx(1:end-1);
%         Sy_end1 = Sy(1:end-1);
%         Sx_end2 = Sx(2:end);
%         Sy_end2 = Sy(2:end); 
%     % Calculate segment lengths and cumulative fault length
%         seglen = sqrt((Sx_end1-Sx_end2).^2+(Sy_end1-Sy_end2).^2);
%         tot_proflen = cumsum(seglen); 
%         add_proflen = [0 tot_proflen'];
%         
        
        
    % end points of general fault
%     xx1 = FData(FAULT).Xfault(1);
%     yy1 = FData(FAULT).Yfault(1);
%     xx2 = FData(FAULT).Xfault(2);
%     yy2 = FData(FAULT).Yfault(2);
        
%% ************************************************************************
%     COMPUTE DISTANCES FROM SURFACE DATA POINTS TO FAULT SEGMENTS
%         ALSO COMPUTE FAULT PARALLEL AND PERPENDICULAR MOTION
%*************************************************************************
% project surface data calculations onto fault segment (FAULT)
% find which segment the surface data points are closest to.
% closest_sec is a vector containing the indeces of the fault segment
% that correspond to each surface datapoint.

% calculate distance to closest segment along fault trace (this
% will capture "shadow zones" at fault bends)
size_grid = size(x_grid);
len_grid = size_grid(1)*size_grid(2);

    closest_seg_data = zeros(size_grid); %zeros(1,len_grid);
    min_dist = closest_seg_data;
    ctrORend = closest_seg_data;
    distalongseg = closest_seg_data;
    westeast = closest_seg_data;
    parslip_grid = closest_seg_data;
    perpslip_grid = closest_seg_data;
    X_project = closest_seg_data;
    Y_project = closest_seg_data;

tic
for n = 1:size_grid(1) % rows
    for m = 1:size_grid(2) % columns
        % project surface data points to fault and store disance and E-W
        X_point = x_grid(n,m);
        Y_point = y_grid(n,m);
        [Min_dist,seg_ind,CtrORend,dist_alongseg,Westeast,x_project,y_project] = project_Kpointstofault_westeast(Sx_end1,Sy_end1,Sx_end2,Sy_end2,X_point,Y_point);
            closest_seg_data(n,m) = seg_ind;
            min_dist(n,m) = Min_dist;
            ctrORend(n,m) = CtrORend;
            distalongseg(n,m) = dist_alongseg;
            westeast(n,m) = Westeast;
            X_project(n,m) = x_project;
            Y_project(n,m) = y_project;
            
        % project NS/ EW slip into fault parallel and perpendicular motion
        if strcmp(DATA_MODE,'global') == 1; 
            az_f = global_fault_az;
            [slip_par,slip_perp,mag_slip] = transform_COSI_faultparperp(az_f,NS_grid(n,m),EW_grid(n,m));
        elseif strcmp(DATA_MODE,'local') == 1; 
            seg_ind = closest_seg_data(n,m); % find index of closest segment
            if seg_ind == 0 % if there is no closest segment use global azimuth
                az_f = global_fault_az;
            else
                az_f = local_fault_az(seg_ind); % otherwise use azimuth of closest segment
            end
            [slip_par,slip_perp,mag_slip] = transform_COSI_faultparperp(az_f,NS_grid(n,m),EW_grid(n,m));
        end
        parslip_grid(n,m) = slip_par;
        perpslip_grid(n,m) = slip_perp;
        % magslip_grid(n,m) = mag_slip;
    end
end
toc

% SURFACE DATA DISTANCE ALONG FAULT
    Surf_distalongfault = distalongseg + add_proflen(closest_seg_data);

% create binary vectors of points east and west of fault
    Surf_east = westeast == 2;
    Surf_west = westeast == 1;

%% ************************************************************************
%               BIN SURFACE DATA BY DISTANCE FROM FAULT
%**************************************************************************
if strcmp(BIN_DIST,'on') == 1;
    % automated bins
    binstart = (0:bin_size:(bin_size*(num_bins-1)));
    binend = (bin_size:bin_size:(bin_size*num_bins));
    zonenum_auto = zeros(size(x_grid));
        for j = 1:num_bins
            binary_find = min_dist >= binstart(j) & min_dist < binend(j);
            zonenum_auto(binary_find) = j;
            zone_distalongfault = Surf_distalongfault(binary_find);
            zone_data_surface = data_surface(binary_find); 
            zone_east = Surf_east(binary_find); % data points on east side of fault that fall within distance range
            zone_west = Surf_west(binary_find);
            
            % collect data for east and west sides of fault
            distalong_west = zone_distalongfault(zone_west);
            distalong_east = zone_distalongfault(zone_east);
            data_surface_west = zone_data_surface(zone_west);
            data_surface_east = zone_data_surface(zone_east); 
            
            % sort data for interpolation (x values must be increasing)
            [sortW,ind_sortW] = sort(distalong_west);
            distalong_w = sortW;
            data_surface_W = data_surface_west(ind_sortW);
            [sortE,ind_sortE] = sort(distalong_east);
            distalong_e = sortE;
            data_surface_E = data_surface_east(ind_sortE);
            
            %% add small number to repeated values for interp (x-values must be distinct)
            % find repeated values and add insignificantly small amount
            % until all values are independent
            sum_find_zero_W = 1;
            distalong_W = distalong_w;
            while sum_find_zero_W ~= 0
                find_same_W = distalong_W(2:end) - distalong_W(1:end-1);
                find_zero_W = find_same_W == 0;
                sum_find_zero_W = sum(find_zero_W);
                if sum_find_zero_W > 0
                    find_zero_add_W = [0 find_zero_W']';
                    distalong_W = distalong_W + (find_zero_add_W)*0.0000001;
                end
            end

            sum_find_zero_E = 1;
            distalong_E = distalong_e;
            while sum_find_zero_E ~= 0
                find_same_E = distalong_E(2:end) - distalong_E(1:end-1);
                find_zero_E = find_same_E == 0;
                sum_find_zero_E = sum(find_zero_E);
                if sum_find_zero_E > 0
                    find_zero_add_E = [0 find_zero_E']';
                    distalong_E = distalong_E + (find_zero_add_E)*0.0000001;
                end
            end
            
            % create smoothed lines through data point clouds
            interp_vect = (1:smooth_win:max(zone_distalongfault))';
            smooth_west = smooth(distalong_W,data_surface_W,smooth_win);
            smooth_east = smooth(distalong_E,data_surface_E,smooth_win);
            smooth_west_interp = interp1(distalong_W,smooth_west,interp_vect);
            smooth_east_interp = interp1(distalong_E,smooth_east,interp_vect);
            net_slip = smooth_west_interp - smooth_east_interp;

            % store values
            Zones(j).binary_find = binary_find;
            Zones(j).G_distalong = zone_distalongfault;
            Zones(j).data_surface = data_surface(binary_find);
            Zones(j).mindist = min_dist(binary_find);
            Zones(j).interp_vect = interp_vect;
            Zones(j).net_slip = net_slip;
            Zones(j).G_west = zone_west; 
            Zones(j).G_east = zone_east;
            Zones(j).smooth_west = smooth_west_interp; 
            Zones(j).smooth_east = smooth_east_interp;
            
        end
end

%% ************************************************************************        
%             DEFINE MAIN DATASET TO USE AS BASE SURFACE
%**************************************************************************
data_surface = parslip_grid;

%% ************************************************************************
%           CALCULATE SLOPES AND CURVATURES FOR SURFACE DATA
%********* *****************************************************************
[slope,slope_az,K,H,P1,P2] = calc_SlopeAzCurvature(data_surface,gridsize);
    % K = Gaussian curvature
    % H =  Mean curvature
    % P1,P2 = Principal curvatures
    
%% ************************************************************************        
%             SAVE DATA FOR QUICK RESTART IN CASE OF CRASH
%**************************************************************************
save Kumamoto_40width_genfault3segs
    
%% ************************************************************************
%   GENERATE PROFILES OF GLENNIE DATA THROUGH FLETCHER FIELD MEASUREMENTS
%**************************************************************************
% cut profiles through each Fletcher data point perpendicular to the
% closest Fletcher segment.   
% prof_len = 1500; % half length of profiles
% prof_width = 20; % half width of profiles 

% RUN PROFILE ANALYZER ONLY IF DIRECTED TO DO SO
if strcmp(PROF_ANALYZE,'on') == 1;
    
    GenerateProfs_SkewNormal_ProfMode_Even
    % OUTPUTS Structure array "Profiles"
        % Profiles.keep_prof - indeces of points in profile
        % Profiles.ctr_easting - profile center points
        % Profiles.ctr_northing
        % Profiles.endpts_x - profile endpoints (x)
        % Profiles.endpts_y
        % Profiles.scarp_alongprof - location of scarps as function of distance along profile
        % Profiles.azimuth_prof -profile azimuth
        % Profiles.dist_alongprof - distance along profile of each surface data point
        % Profiles.accept_results - does user want to keep results for this profile?
        % Profiles.linear_extrap_disp
        % Profiles.linear_extrap_error
        % Profiles.linear_disp_cum
        % Profiles.linear_err_cum
        % Profiles.fit_quality = NaN;
        % RESULTS AUTOMATICALLY SAVED to the folder 'C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Results\Profile_Results\'    

    %% COMPILE ALL SLIP DATA INTO SINGLE STRUCTURE ARRAY
    % if profile offset data is disctibuted among files (due to
    % interruption in code

%     if FAULT == 10 
%         load('C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Results\Profiles_Fault10\Prof_outputs_Mode_even_20170326T200622.mat')
%         Prof1 = Profiles;
%         load('C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Results\Profiles_Fault10\Prof_outputs_Mode_even_20170326T212332.mat')
%         Prof2 = Profiles;
%         load('C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Results\Profiles_Fault10\Prof_outputs_Mode_even_20170326T223409.mat')
%         Prof3 = Profiles;
%         load('C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Results\Profiles_Fault10\Prof_outputs_Mode_even_20170327T150518.mat')
%         Prof4 = Profiles;
%         load('C:\Users\llajoie\Documents\CSM_PhD\PROJECTS\Japan\Results\Profiles_Fault10\Prof_outputs_Mode_even_20170327T233154.mat')
%         Prof5 = Profiles;
% 
%         ProfsALL(40:65) = Prof1(40:65);
%         ProfsALL(66:111) = Prof2(66:111);
%         ProfsALL(112:175) = Prof3(112:175);
%         ProfsALL(176:232) = Prof4(176:232);
%         ProfsALL(233:405) = Prof5(233:405);   
% 
%         distalong_transect = [0:prof_space:(prof_space*(max(size(ProfsALL))-1))];
% 
%             ind_disp = [];
%             ind_error = [];
%             distalong_transect_ind = [];
%             for jj = 40:405
%                 IND_disp = [ProfsALL(jj).linear_extrap_disp];
%                 IND_err = [ProfsALL(jj).linear_extrap_error];
%                     find_Inferr = isinf(IND_err); % remove inf values of error 
%                     IND_Error = IND_err;
%                     IND_Error(find_Inferr) = NaN;
%                 IND_dist = distalong_transect(jj);
%                     IND_distance = ones(size(IND_disp))*IND_dist;
% 
%                 ind_disp = [ind_disp IND_disp'];
%                 ind_error = [ind_error IND_Error'];
%                 distalong_transect_ind = [distalong_transect_ind IND_distance'];
%             end
% 
%             Profs = ProfsALL;
%     end

        %% CREATE SLIP AND DISTANCE VECTORS

        % for error and slip - cumulative (for each data point)
        cum_err = [Profs.linear_err_cum];
        cum_disp = [Profs.linear_disp_cum];
        find_inferr = isinf(cum_err);
        cum_error = cum_err;
        cum_error(find_inferr) = NaN;

        %% PROJECT FIELD MEASUREMENTS ONTO TRANSECT
        X_line = [Profs(end).ctr_easting Profs(PROF_START).ctr_easting];  % or ProfsALL
        Y_line = [Profs(end).ctr_northing Profs(PROF_START).ctr_northing];
        trans_len = sqrt((X_line(1) - X_line(end))^2 + (Y_line(1) - Y_line(end))^2);
        X_points = f_easting;
        Y_points = f_northing;
        [field_xvals,field_yvals,posnegfield,dist_field_alongfault] = project_pointstoline_Kum(X_line,Y_line,X_points,Y_points);
    %         [Min_dist_field,seg_ind_field,CtrORend_field,dist_alongtrans_field,Westeast_field] = project_Kpointstofault_westeast(Sx_end1,Sy_end1,Sx_end2,Sy_end2,X_point,Y_point);
    % end
end 

%% ************************************************************************
%                               FIGURES
%**************************************************************************

%% FIGURE 1 - PLOT GLENNIE DATA FIELD WITH SCARP AND FIELD MEASUREMENTS - MAP
Figure1_DataPlot

%% ONLY RUN THESE FIGURES FOR DISTANCE BIN MODE
if strcmp(BIN_DIST,'on') == 1;
    distbin_plot = 2;
    zone_W = Zones(distbin_plot).G_west; % RED
    zone_E = Zones(distbin_plot).G_east; % BLUE

    %% FIGURE 2 - PLOT HANGING WALL AND FOOTWALL VERTICAL DISPLACEMENTS ALONG CURVING FAULT
    % Figure2_HWFWdisplacementsalongfault
    
    %% FIGURE 3 - PLOT POINTS AS FUNCTION OF DISTANCE FROM FAULT 
    % Figure3_allpointsfunctionofdist
    
    %% FIGURE 4 - PLOT BEST-FIT LINES THROUGH EACH BIN
    % Figure4_bestfitdistbins_kum
    
    %% FIGURE 5 - PLOT POINTS AS FUNCTION OF DISTANCE FROM FAULT - subplots
    % FIGURE 6 - PLOT POINTS IN SELECTED DISTANCE WINDOW AS FUNCTION OF DIST FROM FAULT
    % Figures5and6_distbin_pointsasfunctionofdist

    %% FIGURE 7 - CONFIRM DISTANCE WINDOWING WORKS       
    % FIGURE 8 - CONFIRM LITHOLOGY PROPERLY KEYED TERRAN FLETCHER
    % FIGURE 9 - CHECK EAST-WEST
    % Figures7_9_check_DistWindows_Lith_EW_KUM
end

%% FIGURE 10 - SLOPE MAP
% Figure10_slopemap

%% FIGURES 11 to 13 - CURVATURE MAPS
% Figures11_13_curvaturemaps
   
%% ONLY RUN THESE FIGURES IF PROFILE ANALYZER IS TURNED ON
if strcmp(PROF_ANALYZE,'on') == 1;
    %% FIGURE 17 - HISTOGRAM OF GAUSS FIT B-VALUES
%     Figure17_bgaussHistogram

    %% FIGURE 18 - RUPTURE ZONE ASYMMETRY
    % FIGURE 20 - RUPTURE ZONE WIDTH
    % Figures18and20_RZasymmetry_RZwidth

    %% FIGURE 19 - DIP VAL v B_GAUSS and RZ WIDTH v C_GAUSS
    % Figure19_DipValBGauss_RZwidthCGauss
    
    %% FIGURE 21 - Profiles - embedded in GenerateProfs_GlennieFletcher_SkewNormal

    %% FIGURES 22-25 - PROFILE STATISTICS
    Figures22_26_ProfileStatsResults_Kum
    
end

%% FIGURE 20 - DISPLACEMENT VECTORS
% Figure20_DisplacementVectors



%         plot_northing = [ProfsALL(5:55).ctr_northing]
%         plot_easting = [ProfsALL(5:55).ctr_easting]

    