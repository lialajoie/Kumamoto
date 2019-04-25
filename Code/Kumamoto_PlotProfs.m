% Lia Lajoie
% CSM
% 2 Nov 2018
% plot Kumamoto profile locations for 3-segment generalized fault

FAULT = 11; % three-segment generalized fault
gridsize = 8;
prof_len = 2500; % half-length profile
prof_width = 20; % half-width of profiles
prof_space = prof_width*2;
DATA_MODE = 'local';
global_fault_az = 45;

% load data
FaultData_Kumamoto_store;
loadfile = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/COSIdata_fault11.mat';
    load(loadfile);
        % NS_grid/ EW_grid: NS/ EW displacements computed by COSI-Corr
        % SNR_grid: signal to noise ratio 
        % x_grid, y_grid, x_vect, y_vect
        
%% ************************************************************************
%                   CALCULATE FAULT SEGMENT DATA
%**************************************************************************      
% create vectors of fault endpoints (x and y values along fault trace)
    Sxx = FaultData(FAULT).Xfault;
    Syy = FaultData(FAULT).Yfault;
    numsegs_FAULT = length(Sxx)-1; % number of fault segments

% ORDER POINTS NORTH TO SOUTH - ALL CALCULATIONS ASSUME THIS POINT ORDER
    [Sy,ind] = sort(Syy,'descend');
    Sx = Sxx(ind);
   
% make vector of start and end points for FAULT segments
    Sx_end1 = Sx(1:end-1);
    Sy_end1 = Sy(1:end-1);
    Sx_end2 = Sx(2:end);
    Sy_end2 = Sy(2:end);
  
% calculate dx, dy, scale factor (a), and length for each fault segment
    dx = Sx_end2-Sx_end1; % dx - vector with length = number of segments
    dy = Sy_end2-Sy_end1; % dy
    a_seg = (prof_space^2./(dx.^2 + dy.^2)).^(1/2); % reduction factor to convert the vector dxdy to a length of profile_separation
    seglen = sqrt((dx).^2+(dy).^2); % length of each segment
    m = dy./dx; % slope of each segment
    norm_m = -(m).^-1; % normal to each segment
    
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
            normseg(profs) = norm_m(jj); % normal of each profile in segment
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
        Px = zeros(numprofs,1);
        Py = Px;
        for kk = 1:numprofs
            seg = profsegs(kk);
            Px(kk) = Sx(seg) + a_start(seg)*dx(seg) + (kk-startprof(seg))*a_seg(seg)*dx(seg); % x-coordinate of profile center
            Py(kk) = Sy(seg) + a_start(seg)*dy(seg) + (kk-startprof(seg))*a_seg(seg)*dy(seg); % y-coordinate of profile center
        end
        
      % Profile endpoints
      dx_prof = ones(size(normseg));
      dy_prof = normseg;
      a_prof = (prof_len^2./(dx_prof.^2 + dy_prof.^2)).^(1/2);
      
        P_end1_x = Px + dx_prof.*a_prof;
        P_end2_x = Px - dx_prof.*a_prof;
        P_end1_y = Py + dy_prof.*a_prof;
        P_end2_y = Py - dy_prof.*a_prof;

%     % Calculate segment lengths and cumulative fault length
        seglen = sqrt((Sx_end1-Sx_end2).^2+(Sy_end1-Sy_end2).^2);
        tot_proflen = cumsum(seglen); 
        add_proflen = [0 tot_proflen];

%% ************************************************************************
%                          IMPORT FIELD DATA
%**************************************************************************
% SPECIFIC DATA POINTS TO FAULT (JUST X,Y)
file_name = '/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/FieldData_Fault11.csv';
    field_fault_data = importdata(file_name);
    x_field_fault = field_fault_data(:,1);
    y_field_fault = field_fault_data(:,2);

% % FIELD DATA - ENTIRE DATASET, NOT FAULT SPECIFIC
%     ParseField_completedataset;
%     % Field_Rslip; Field_Sup; Field_dipslip; Field_northing; Field_easting
%    
% % CROSS REFERENCE FIELD DATA
% % Cross reference field data. I  plotted all field measurements
% % in Matlab and filtered datapoints to only those located on each specific 
% % fault trace. When saving the data, only x and y locations are retained, 
% % so this step cross-references the x,y locations of field data on each 
% % fault to the full dataset and extracts other data values.
%     field_indeces = [];
%     for t = 1:length(x_field_fault)
%         x_find = find(Field_easting == x_field_fault(t) & Field_northing == y_field_fault(t));
%         field_indeces = [field_indeces x_find'];
%     end
%         
% % COLLECT FAULT SPECIFIC FIELD DATA
%     f_easting = Field_easting(field_indeces); 
%     f_northing = Field_northing(field_indeces); 
%     f_Rslip = Field_Rslip(field_indeces); 
%     f_Sup = Field_Sup(field_indeces); 
%     f_dipslip = Field_dipslip(field_indeces); 
        
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
    
%% FIGURE 1 
        figure(1)
        clf
        for ii = 1:5:length(Px)
            plot([P_end1_x(ii) P_end2_x(ii)],[P_end1_y(ii) P_end2_y(ii)],'k')
            hold on
        end
        plot(Px,Py,'ro');
        hold on
        plot(x_field_fault,y_field_fault,'ko')

        axis equal

%% FIGURE 2

Sx_plot = [-9257 -11760 -18600 -19576];
Sy_plot = [-19459 -20496 -27000 -30264];
D_plot = ((Sx_plot(2:end) - Sx_plot(1:end-1)).^2 + (Sy_plot(2:end) - Sy_plot(1:end-1)).^2).^(1/2);

c_range = [-2 2];
[blueredramp] = make_blueredramp(c_range,'zero');
figure(2)
clf
   s = surf(x_grid,y_grid,parslip_grid);
   hold on
    for ii = 20:5:408
        plot3([P_end1_x(ii) P_end2_x(ii)],[P_end1_y(ii) P_end2_y(ii)],[100 100],'k','linewidth',2)
        hold on
    end
    hold on
    plot3(Sx_plot,Sy_plot,ones(size(Sy))*100,'k-','linewidth',3)
       c = colorbar
       c.Label.String = 'pixel displacement in direction parallel to closest transect segment (m)';
       s.EdgeColor = 'none';
       colormap(blueredramp)
       caxis([-2 2])
       ax.XRuler.Exponent = 0;
       ax.YRuler.Exponent = 0;
       ax.XTickLabelRotation = 45;
       ax.YTickLabelRotation = 45;

%        xlim([-2.45e04 -0.65e04])
       view(2)
       axis equal
       grid off
       xlabel('easting (m)')
       ylabel('northing (m)')
       set(gca,'fontsize',16,'fontweight','bold')
       
       
       
       
       