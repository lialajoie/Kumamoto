% Project field measurements onto fault

Field = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/FieldData/KumFieldData_plot_UPDATED.mat');
fx = Field.f_easting;
fy = Field.f_northing;
slip = Field.f_Rslip

tx = [-19394 -18603 -11768 -9257];
ty = [-29656 -27011 -20503 -19459];



X_line1 = tx(1:end-1); X_line2 = tx(2:end);
Y_line1 = ty(1:end-1); Y_line2 = ty(2:end);
seglen = sqrt((X_line1-X_line2).^2 + (Y_line1-Y_line2).^2);
seg_add = [0 cumsum(seglen(1:end-1))];
dist_along_seg = zeros(size(fx));
dist_along_fault = zeros(size(fx));
seg_indeces = zeros(size(fx));
for ii = 1:length(fx)
    X_point = fx(ii);
    Y_point = fy(ii);
    [min_dist,seg_ind,ctrORend,dist_alongseg,westeast] = project_pointstofault_westeast_Kumamoto(X_line1,Y_line1,X_line2,Y_line2,X_point,Y_point)
    if ctrORend == 1
        d_seg = seglen(seg_ind) - dist_alongseg
        dist_along_seg(ii) = d_seg;
        dist_along_fault(ii) = d_seg + seg_add(seg_ind);
        seg_indeces(ii) = seg_ind;
    end
end


figure(1)
clf
for kk = 1:length(fx)
    if seg_indeces(kk) == 1
        scatter(fx(kk),fy(kk),'r','filled')
    elseif seg_indeces(kk) == 2
        scatter(fx(kk),fy(kk),'g','filled')
    elseif seg_indeces(kk) == 3
        scatter(fx(kk),fy(kk),'b','filled')
    end 
    hold on
end
hold on
plot(tx,ty,'k-')
colormap(jet)


figure(3)
clf
for kk = 1:length(fx)
    if seg_indeces(kk) == 1
        scatter(dist_along_fault(kk),slip(kk),'r','filled')
    elseif seg_indeces(kk) == 2
        scatter(dist_along_fault(kk),slip(kk),'g','filled')
    elseif seg_indeces(kk) == 3
        scatter(dist_along_fault(kk),slip(kk),'b','filled')
    end 
    hold on
end
scatter(Field.F_distalongfault,Field.f_Rslip,'k')
hold on
plot([seg_add(2) seg_add(2)],[-100 250],'k--','linewidth',2)
hold on
plot([seg_add(3) seg_add(3)],[-100 250],'k--','linewidth',2)

Field.F_distalongfault = [];
Field.F_distalongfault = dist_along_fault;
figure(4)
scatter(Field.F_distalongfault,Field.f_Rslip)

save('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/FieldData/KumFieldData_plot_FINAL.mat','Field')
