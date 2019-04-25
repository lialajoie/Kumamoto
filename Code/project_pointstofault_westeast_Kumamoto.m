function [min_dist,seg_ind,ctrORend,dist_alongseg,westeast] = project_pointstofault_westeast_Kumamoto(X_line1,Y_line1,X_line2,Y_line2,X_point,Y_point)
    % compute distance from point to each fault segment using: 
    % distance = (a*x0 + b*y0 + c)/sqrt(a^2 + b^2) for 0 = ax + by + c.
    % When converting from y = mx + b, a = m = sslope, b = 1, c = b = bval.
%     X_point = x_grid(100,100);
%     Y_point = y_grid(100,100);
%         
%     X_line1 = Sx_end1;
%     Y_line1 = Sy_end1;
%     X_line2 = Sx_end2;
%     Y_line2 = Sy_end2;
    
    Min_dist = zeros(size(X_line1));
    CtrORend = zeros(size(X_line1));
    Dist_pt2norm = zeros(size(X_line1));
    Dist_alongseg = zeros(size(X_line1));
    WestEast = zeros(size(X_line1));
    for i = 1:length(X_line1) % for each line segment
        % end points of line segment (i)
            x1 = X_line1(i);
            y1 = Y_line1(i);
            x2 = X_line2(i);
            y2 = Y_line2(i);
        
        % Only execute for segments with non-zero length
        if x1 == x2 && y1 == y2
            Min_dist(i) = NaN;
            Dist_pt2norm(i) = NaN;
            Dist_alongseg(i) = NaN;
            CtrORend(i) = NaN;
        else
            % always order points such that Y1 is greater than Y2
                if y1 > y2
                    X1 = x1;
                    Y1 = y1;
                    X2 = x2;
                    Y2 = y2;
                else
                    X1 = x2;
                    Y1 = y2;
                    X2 = x1;
                    Y2 = y1;
                end

            % calculate slope and b-value (for y = mx + b) of segment (i)
                dx = X2-X1;
                dy = Y2-Y1;
                slope = dy/dx;
                b = Y1 - slope*X1;

            % calculate distance from point to line segment (i) and
            % determine whether point lies east or west of fault
                dist_P2seg = (-slope*X_point + Y_point - b)/(sqrt(slope^2+1));
                abs_dist_P2seg = abs(dist_P2seg);
                
                % West = 1, East = 2
%                 if slope >= 0 && slope < 1 % for shallow, positive slopes, positive distance is east
%                     find_EW = (dist_P2seg >= 0);
%                     if find_EW > 0
%                         WestEast(i) = 2;
%                     else
%                         WestEast(i) = 1;
%                     end
                if slope >= 0 % for positive slopes, positive distance is west
                    find_EW = (dist_P2seg >= 0);
                    if find_EW > 0
                        WestEast(i) = 1;
                    else
                        WestEast(i) = 2;
                    end
                elseif slope < 0 % for negative slopes, positive distance is east
                    find_EW = (dist_P2seg >= 0);
                    if find_EW > 0
                        WestEast(i) = 2;
                    else
                        WestEast(i) = 1;
                    end
                end

            % project data point orthogonally onto segment and calculate distance
            % from projected point to start of segment (for later use in
            % distance along fault/ profile)
                a = sqrt(dist_P2seg.^2/(dx^2 + dy^2));
                    pos_ind = dist_P2seg <= 0;
                    neg_ind = dist_P2seg > 0;
                posneg = pos_ind + neg_ind*(-1);
                A = a.*posneg;
                x_project = X_point - A*dy;
                y_project = Y_point + A*dx;
                

            %% determine whether points are contained within normal space to segment
            % keep_pts is logical with same length as condfind
            % keep_indeces converts logical to indeces referenced to original dataset
            % execute by calculating distance of each point in segment from lines
            % running perpendicular to segment ends, use positive and negative
            % distances to determine which are within segment.
                inv_slope = -dx/dy;
                inv_bval_1 = Y1 - inv_slope*X1;
                inv_bval_2 = Y2 - inv_slope*X2;
                dist_data2L_1 = (-inv_slope*X_point + Y_point - inv_bval_1)/(sqrt(inv_slope^2+1));
                dist_data2L_2 = (-inv_slope*X_point + Y_point - inv_bval_2)/(sqrt(inv_slope^2+1));

            % decide if data point lies within a segment, or outside of it
            % Teran end point 1 is always north of end point 2 in my computations
                within_seg = dist_data2L_1 <= 0 & dist_data2L_2 >= 0;

            if within_seg == 1 % if point is normal to a fault segment
                Min_dist(i) = abs_dist_P2seg;
                CtrORend(i) = 1;
                Dist_pt2norm(i) = 0;
                Dist_alongseg(i) = sqrt((X1-x_project)^2 + (Y1-y_project)^2);
            else % if point falls in "dead space" on outer fault bend
                % find distance to segment endpoints (Preserving original point
                % order)
                    d1 = sqrt((x1-X_point)^2+(y1-Y_point)^2);
                    d2 = sqrt((x2-X_point)^2+(y2-Y_point)^2);
                % store distance to closer end point and to normal line at endpoint 
                % used to bisect angle between normal lines of adjacent segments
                if d1 < d2 % if the point is closer to endpoint 1
                    Min_dist(i) = d1;
                    Dist_pt2norm(i) = abs(dist_data2L_1);
                    Dist_alongseg(i) = 0;
                elseif d2 < d1 % if the point is closet to endpoint 
                    Min_dist(i) = d2;
                    Dist_pt2norm(i) = abs(dist_data2L_2);
                    Dist_alongseg(i) = sqrt((X1-X2)^2 + (Y1-Y2)^2);
                end
                CtrORend(i) = 2;
            end
        end
    end
    
    %% find closest line segment
    [sort_val,sort_ind] = sort(Min_dist);
    ctrORend_sort1 = CtrORend(sort_ind(1));
    
    % if the closest distance to the fault is at a segment end
    if ctrORend_sort1 == 2
        dn1 = Dist_pt2norm(sort_ind(1));
        dn2 = Dist_pt2norm(sort_ind(2));
        if dn1 <= dn2
            min_dist = Min_dist(sort_ind(1));
            seg_ind = sort_ind(1);
            dist_alongseg = Dist_alongseg(sort_ind(1));
            ctrORend = CtrORend(sort_ind(1));
            westeast = WestEast(sort_ind(1));
        elseif dn2 < dn1
            min_dist = Min_dist(sort_ind(2));
            seg_ind = sort_ind(2);
            dist_alongseg = Dist_alongseg(sort_ind(2));
            ctrORend = CtrORend(sort_ind(2));
            westeast = WestEast(sort_ind(2));
        end
    % if the closest distance to the fault is at a segment center
    else
        min_dist = sort_val(1);
        seg_ind = sort_ind(1);
        dist_alongseg = Dist_alongseg(sort_ind(1));
        ctrORend = CtrORend(sort_ind(1));
        westeast = WestEast(sort_ind(1));
    end

end

