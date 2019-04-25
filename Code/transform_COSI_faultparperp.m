function [slip_par,slip_perp,mag_slip] = transform_COSI_faultparperp(az_f,NS,EW)
    % az_f = fault azimuth
    % az_n = fault_normal direction
    % NS/ EW = NS/ EW cosi corr displacements
    az_n = 360-(90-az_f); % normal azimuth
    %% CALCULATE SLIP PARAMETERS
% 1) caulculate angle of slip with respect to "horizontal" or due east
    theta = abs(atand(NS./EW)); % angle of slip vector from horizontal
% 2) determine which quadrant the vector is in (upper right is pos/pos)
    Q1 = NS >= 0 & EW >= 0; % quadrant 1 - upper right
    Q2 = NS >= 0 & EW < 0; % upper left
    Q3 = NS < 0 & EW < 0; % lower left
    Q4 = NS < 0 & EW >= 0; % lower right
% 3) calculate azimuth of slip
    slip_az = (90-theta).*Q1;
    slip_az = slip_az + (270+theta).*Q2;
    slip_az = slip_az + (270-theta).*Q3;
    slip_az = slip_az + (90+theta).*Q4;
% 4) calculate magnitude of slip
    mag_slip = sqrt(NS.^2 + EW.^2);
    
%% TRANSFORM TO ALONG FAULT DIRECTION
% 5) determine angle relative to azimuth
    diff_par = abs(slip_az-az_f);
% 6) project slip onto azimuth
    slip_par = mag_slip.*cosd(diff_par);
    
%% TRANSFORM TO FAULT PERPENDICULAR DIRECTION
% 5) determine angle relative to azimuth
    diff_perp = abs(slip_az-az_n);
% 6) project slip onto azimuth
    slip_perp = mag_slip.*cosd(diff_perp);
    
end