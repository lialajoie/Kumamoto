% Lia Lajoie
% CSM 27 Oct 2018
% plot all RMS v window length data for all EQ on one plot
clear all

B = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Balochistan_plot_data.mat');
    B_win_size = B.win_size; B_RMS_win = B.RMS_win; B_distalong_transect_m = B.distalong_transect_m;
    B_sigma_RMS_dt = B.sigma_RMS_dt; B_RMS_win_sigma = B.RMS_win_sigma; B_H_RMSslope_allpoints = B.H_RMSslope_allpoints;
    B_H_RMSslope_limited = B.H_RMSslope_limited; B_H_RMSslope_limited2 = B.H_RMSslope_limited2;
Bz = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Balu_Zinke_plot_data.mat');
    Bz_win_size = Bz.win_size; Bz_RMS_win = Bz.RMS_win; Bz_distalong_transect_m = Bz.distalong_transect_m;
    Bz_sigma_RMS_dt = Bz.sigma_RMS_dt; Bz_RMS_win_sigma = Bz.RMS_win_sigma; Bz_H_RMSslope_allpoints = Bz.H_RMSslope_allpoints;
    Bz_H_RMSslope_limited = Bz.H_RMSslope_limited; Bz_H_RMSslope_limited2 = Bz.H_RMSslope_limited2;
HM = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/HectorMine_plot_data.mat');
    HM_win_size = HM.win_size; HM_RMS_win = HM.RMS_win; HM_distalong_transect_m = HM.distalong_transect_m;
    HM_sigma_RMS_dt = HM.sigma_RMS_dt; HM_RMS_win_sigma = HM.RMS_win_sigma; HM_H_RMSslope_allpoints = HM.H_RMSslope_allpoints;
    HM_H_RMSslope_limited = HM.H_RMSslope_limited; HM_H_RMSslope_limited2 = HM.H_RMSslope_limited2;
K_all = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Kumamoto_plot_data_FAULT8.mat');
    Ka_win_size = K_all.win_size; Ka_RMS_win = K_all.RMS_win; Ka_distalong_transect_m = K_all.DIST;
    Ka_sigma_RMS_dt = K_all.sigma_RMS_dt; Ka_RMS_win_sigma = K_all.RMS_win_sigma;
    Ka_H_RMSslope_allpoints = K_all.H_RMSslope_allpoints; Ka_H_RMSslope_limited = K_all.H_RMSslope_limited;
K1 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Kumamoto_plot_data_FAULT1.mat');
    K1_win_size = K1.win_size; K1_RMS_win = K1.RMS_win; K1_distalong_transect_m = K1.DIST;
    K1_sigma_RMS_dt = K1.sigma_RMS_dt; K1_RMS_win_sigma = K1.RMS_win_sigma;
    K1_H_RMSslope_allpoints = K1.H_RMSslope_allpoints; K1_H_RMSslope_limited = K1.H_RMSslope_limited;
K2 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Kumamoto_plot_data_FAULT2.mat');
    K2_win_size = K2.win_size; K2_RMS_win = K2.RMS_win; K2_distalong_transect_m = K2.DIST;
    K2_sigma_RMS_dt = K2.sigma_RMS_dt; K2_RMS_win_sigma = K2.RMS_win_sigma;
    K2_H_RMSslope_allpoints = K2.H_RMSslope_allpoints; K2_H_RMSslope_limited = K2.H_RMSslope_limited;
K4 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Kumamoto_plot_data_FAULT4.mat');
    K4_win_size = K4.win_size; K4_RMS_win = K4.RMS_win; K4_distalong_transect_m = K4.DIST;
    K4_sigma_RMS_dt = K4.sigma_RMS_dt; K4_RMS_win_sigma = K4.RMS_win_sigma;
    K4_H_RMSslope_allpoints = K4.H_RMSslope_allpoints; K4_H_RMSslope_limited = K4.H_RMSslope_limited;
K5 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Kumamoto_plot_data_FAULT5.mat');
    K5_win_size = K5.win_size; K5_RMS_win = K5.RMS_win; K5_distalong_transect_m = K5.DIST;
    K5_sigma_RMS_dt = K5.sigma_RMS_dt; K5_RMS_win_sigma = K5.RMS_win_sigma;
    K5_H_RMSslope_allpoints = K5.H_RMSslope_allpoints; K5_H_RMSslope_limited = K5.H_RMSslope_limited;
L = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Landers_plot_data.mat');
    L_win_size = L.win_size; L_RMS_win = L.RMS_win; L_distalong_transect_m = L.distalong_transect_m;
    L_sigma_RMS_dt = L.sigma_RMS_dt; L_RMS_win_sigma = L.RMS_win_sigma;
    L_H_RMSslope_allpoints = L.H_RMSslope_allpoints; L_H_RMSslope_limited = L.H_RMSslope_limited;
E1 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/EMC_fault1_plot_data.mat');
    E1_win_size = E1.win_size; E1_RMS_win = E1.RMS_win; E1_distalong_transect_m = E1.distalong_transect_m;
    E1_sigma_RMS_dt = E1.sigma_RMS_dt; E1_RMS_win_sigma = E1.RMS_win_sigma;
    E1_H_RMSslope_allpoints = E1.H_RMSslope_allpoints; E1_H_RMSslope_limited = E1.H_RMSslope_limited;
E2 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/EMC_fault2_plot_data.mat');
    E2_win_size = E2.win_size; E2_RMS_win = E2.RMS_win; E2_distalong_transect_m = E2.distalong_transect_m;
    E2_sigma_RMS_dt = E2.sigma_RMS_dt; E2_RMS_win_sigma = E2.RMS_win_sigma;
    E2_H_RMSslope_allpoints = E2.H_RMSslope_allpoints; E2_H_RMSslope_limited = E2.H_RMSslope_limited;
E3 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/EMC_fault3_plot_data.mat');
    E3_win_size = E3.win_size; E3_RMS_win = E3.RMS_win; E3_distalong_transect_m = E3.distalong_transect_m;
    E3_sigma_RMS_dt = E3.sigma_RMS_dt; E3_RMS_win_sigma = E3.RMS_win_sigma;
    E3_H_RMSslope_allpoints = E3.H_RMSslope_allpoints; E3_H_RMSslope_limited = E3.H_RMSslope_limited;
E45 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/EMC_fault5_plot_data.mat');
    E45_win_size = E45.win_size; E45_RMS_win = E45.RMS_win; E45_distalong_transect_m = E45.distalong_transect_m;
    E45_sigma_RMS_dt = E45.sigma_RMS_dt; E45_RMS_win_sigma = E45.RMS_win_sigma;
    E45_H_RMSslope_allpoints = E45.H_RMSslope_allpoints; E45_H_RMSslope_limited = E45.H_RMSslope_limited;
E8 = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/EMC_DipSlip/EMC_fault8_plot_data.mat');
    E8_win_size = E8.win_size; E8_RMS_win = E8.RMS_win; E8_distalong_transect_m = E8.distalong_transect_m;
    E8_sigma_RMS_dt = E8.sigma_RMS_dt; E8_RMS_win_sigma = E8.RMS_win_sigma;
    E8_H_RMSslope_allpoints = E8.H_RMSslope_allpoints; E8_H_RMSslope_limited = E8.H_RMSslope_limited;
P = open('/Users/llajoie/Documents/CSM_PhD/Kumamoto/Data/SlipProfiles/Palu_plot_data.mat');
    P_win_size = P.win_size; P_RMS_win = P.RMS_win; P_distalong_transect_m = P.distalong_transect_m;
    P_sigma_RMS_dt = P.sigma_RMS_dt; P_RMS_win_sigma = P.RMS_win_sigma;
    P_H_RMSslope_allpoints = P.H_RMSslope_allpoints; P_H_RMSslope_limited = P.H_RMSslope_limited;
    
%% FIGURE 1 - HURST FROM RMS ROUGHNESS, ALL EARTHQUAKES
% background grey
    c10 = [0.8 0.8 0.8];
% kumamoto colors
    cKa = 'b';
    cK1 = [102 51 153]/256;
    cK2 = [255 102 0]/256;
    cK4 = [51 153 0]/256;
    cK5 = [1, 0, 1];
% EMC colors
    cE1 = [102 255 0]/256;
    cE2 = [0 255 255]/256;
    cE3 = [0 0 0]/256;
    cE45 = [255 153 51]/256;
    cE8 = [153 51 255]/256;
% other fault colors
    cB = [255 204 0]/256;
    cBz = [255 102 153]/256;
    cHM = [0 204 255]/256;
    cL = 'r'; 
    cP = [153 204 000]/256;
    
    
    %neww cmment.


figure(1)
clf
    subplot(3,1,1)
    % plot Balochistan
        plot(B_win_size,B_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(B_distalong_transect_m(end)/1000,B_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(B_win_size,B_RMS_win,B_RMS_win_sigma,'color',c10)
%         hold on
    % plot Palu
        plot(P_win_size,P_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(P_distalong_transect_m(end)/1000,P_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
    % plot Hector Mine
        plot(HM_win_size,HM_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(HM_distalong_transect_m(end)/1000,HM_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(HM_win_size,HM_RMS_win,HM_RMS_win_sigma,'color',c10)
%         hold on
    % plot Landers
        plot(L_win_size,L_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(L_distalong_transect_m(end)/1000,L_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(L_win_size,L_RMS_win,L_RMS_win_sigma,'color',c10)
%         hold on
    % plot EMC - P, LS
        plot(E1_win_size,E1_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(E1_distalong_transect_m(end)/1000,E1_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E1_win_size,E1_RMS_win,E1_RMS_win_sigma,'color',c10)
%         hold on
    % plot EMC - Borrego
        plot(E2_win_size,E2_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E2_distalong_transect_m(end)/1000,E2_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E2_win_size,E2_RMS_win,E2_RMS_win_sigma,'color',c10)
%         hold on
    % plot EMC - PI
        plot(E3_win_size,E3_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E3_distalong_transect_m(end)/1000,E3_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
        errorbar(E3_win_size,E3_RMS_win,E3_RMS_win_sigma,'color',c10) 
        hold on
    % plot EMC - PS
        plot(E45_win_size,E45_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E45_distalong_transect_m(end)/1000,E45_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E45_win_size,E45_RMS_win,E45_RMS_win_sigma,'color',c10)   
%         hold on
    % plot EMC - Puerta
        plot(E8_win_size,E8_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E8_distalong_transect_m(end)/1000,E8_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E8_win_size,E8_RMS_win,E8_RMS_win_sigma,'color',c10)
%         hold on
    % plot Kumamoto - all
        plot(Ka_win_size,Ka_RMS_win,'color',cKa,'linewidth',2)
        hold on
        plot((max(Ka_distalong_transect_m)-min(Ka_distalong_transect_m)),Ka_sigma_RMS_dt,'*','color',cKa,'MarkerSize',12)
        hold on
        errorbar(Ka_win_size,Ka_RMS_win,Ka_RMS_win_sigma,'color',cKa)
        hold on
    % plot K1
        plot(K1_win_size,K1_RMS_win,'color',cK1,'linewidth',2)
        hold on
        plot((max(K1_distalong_transect_m)-min(K1_distalong_transect_m)),K1_sigma_RMS_dt,'*','color',cK1,'MarkerSize',12)
        hold on
        errorbar(K1_win_size,K1_RMS_win,K1_RMS_win_sigma,'color',cK1)
        hold on
    % plot K2
        plot(K2_win_size,K2_RMS_win,'color',cK2,'linewidth',2)
        hold on
        plot((max(K2_distalong_transect_m)-min(K2_distalong_transect_m)),K2_sigma_RMS_dt,'*','color',cK2,'MarkerSize',12)
        hold on
        errorbar(K2_win_size,K2_RMS_win,K2_RMS_win_sigma,'color',cK2)
        hold on
    % plot K4
        plot(K4_win_size,K4_RMS_win,'color',cK4,'linewidth',2)
        hold on
        plot((max(K4_distalong_transect_m)-min(K4_distalong_transect_m)),K4_sigma_RMS_dt,'*','color',cK4,'MarkerSize',12)
        hold on
        errorbar(K4_win_size,K4_RMS_win,K4_RMS_win_sigma,'color',cK4)
        hold on
    % plot K5
        plot(K5_win_size,K5_RMS_win,'color',cK5,'linewidth',2)
        hold on
        plot((max(K5_distalong_transect_m)-min(K5_distalong_transect_m)),K5_sigma_RMS_dt,'*','color',cK5,'MarkerSize',12)
        hold on
        errorbar(K5_win_size,K5_RMS_win,K5_RMS_win_sigma,'color',cK5)
        
        
        % make legend
            h1(1) = plot(NaN,NaN,'-','color',cKa,'linewidth',2);
            h1(2) = plot(NaN,NaN,'-','color',cK1,'linewidth',2);
            h1(3) = plot(NaN,NaN,'-','color',cK2,'linewidth',2);
            h1(4) = plot(NaN,NaN,'-','color',cK4,'linewidth',2);
            h1(5) = plot(NaN,NaN,'-','color',cK5,'linewidth',2);
            legend(h1,'All segments','Segment 1','Segment 2',...
                'Segment 4','Segment 5')

            xlim([3*10^-1 3*10^2])
            ylim([7*10^-2 4])
            title('2016 Kumamoto earthquake')
            ax = gca;
            ax.XTickLabel = [];
%             ylabel('average value of RMS variance (m)')
            set(gca,'fontsize',16,'fontweight','bold','yscale','log','xscale','log')
            grid on
            
    subplot(3,1,2)
    % plot Balochistan
        plot(B_win_size,B_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(B_distalong_transect_m(end)/1000,B_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
    % plot Palu
        plot(P_win_size,P_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(P_distalong_transect_m(end)/1000,P_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(B_win_size,B_RMS_win,B_RMS_win_sigma,'color',c10)
%         hold on
    % plot Hector Mine
        plot(HM_win_size,HM_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(HM_distalong_transect_m(end)/1000,HM_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(HM_win_size,HM_RMS_win,HM_RMS_win_sigma,'color',c10)
%         hold on
    % plot Landers
        plot(L_win_size,L_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(L_distalong_transect_m(end)/1000,L_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(L_win_size,L_RMS_win,L_RMS_win_sigma,'color',c10)
%         hold on
    % plot Kumamoto - all
        plot(Ka_win_size,Ka_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(Ka_distalong_transect_m)-min(Ka_distalong_transect_m)),Ka_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(Ka_win_size,Ka_RMS_win,Ka_RMS_win_sigma,'color',c10)
%         hold on
    % plot K1
        plot(K1_win_size,K1_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K1_distalong_transect_m)-min(K1_distalong_transect_m)),K1_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K1_win_size,K1_RMS_win,K1_RMS_win_sigma,'color',c10)
%         hold on
    % plot K2
        plot(K2_win_size,K2_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K2_distalong_transect_m)-min(K2_distalong_transect_m)),K2_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K2_win_size,K2_RMS_win,K2_RMS_win_sigma,'color',c10)
%         hold on
    % plot K4
        plot(K4_win_size,K4_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K4_distalong_transect_m)-min(K4_distalong_transect_m)),K4_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K4_win_size,K4_RMS_win,K4_RMS_win_sigma,'color',c10)
%         hold on
    % plot K5
        plot(K5_win_size,K5_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K5_distalong_transect_m)-min(K5_distalong_transect_m)),K5_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K5_win_size,K5_RMS_win,K5_RMS_win_sigma,'color',c10)
%         hold on
    % plot EMC - P, LS
        plot(E1_win_size,E1_RMS_win,'-','color',cE1,'linewidth',2)
        hold on
        plot(E1_distalong_transect_m(end)/1000,E1_sigma_RMS_dt,'*','color',cE1,'MarkerSize',12)
        hold on
        errorbar(E1_win_size,E1_RMS_win,E1_RMS_win_sigma,'color',cE1)
        hold on
    % plot EMC - Borrego
        plot(E2_win_size,E2_RMS_win,'color',cE2,'linewidth',2)
        hold on
        plot(E2_distalong_transect_m(end)/1000,E2_sigma_RMS_dt,'*','color',cE2,'MarkerSize',12)
        hold on
        errorbar(E2_win_size,E2_RMS_win,E2_RMS_win_sigma,'color',cE2)
        hold on
    % plot EMC - PI
        plot(E3_win_size,E3_RMS_win,'color',cE3,'linewidth',2)
        hold on
        plot(E3_distalong_transect_m(end)/1000,E3_sigma_RMS_dt,'*','color',cE3,'MarkerSize',12)
        hold on
        errorbar(E3_win_size,E3_RMS_win,E3_RMS_win_sigma,'color',cE3) 
        hold on
    % plot EMC - PS
        plot(E45_win_size,E45_RMS_win,'color',cE45,'linewidth',2)
        hold on
        plot(E45_distalong_transect_m(end)/1000,E45_sigma_RMS_dt,'*','color',cE45,'MarkerSize',12)
        hold on
        errorbar(E45_win_size,E45_RMS_win,E45_RMS_win_sigma,'color',cE45)   
        hold on
    % plot EMC - Puerta
        plot(E8_win_size,E8_RMS_win,'color',cE8,'linewidth',2)
        hold on
        plot(E8_distalong_transect_m(end)/1000,E8_sigma_RMS_dt,'*','color',cE8,'MarkerSize',12)
        hold on
        errorbar(E8_win_size,E8_RMS_win,E8_RMS_win_sigma,'color',cE8)
        hold on
        
        
        % make legend
            h2(1) = plot(NaN,NaN,'-','color',cE1,'linewidth',2);
            h2(2) = plot(NaN,NaN,'-','color',cE2,'linewidth',2);
            h2(3) = plot(NaN,NaN,'-','color',cE3,'linewidth',2);
            h2(4) = plot(NaN,NaN,'-','color',cE45,'linewidth',2);
            h2(5) = plot(NaN,NaN,'-','color',cE8,'linewidth',2);
            legend(h2,'Pescadores, Laguna Salada','Borrego','Paso Inferior',...
                'Paso Superior','Puerta accomodation zone')

            xlim([3*10^-1 3*10^2])
            ylim([7*10^-2 4])
            title('2010 El Mayor-Cucapah earthquake')
            ax = gca;
            ax.XTickLabel = [];
            ylabel('average value of RMS variance (m)')
            set(gca,'fontsize',16,'fontweight','bold','yscale','log','xscale','log')
            grid on
        
    subplot(3,1,3)
    % plot Kumamoto - all
        plot(Ka_win_size,Ka_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(Ka_distalong_transect_m)-min(Ka_distalong_transect_m)),Ka_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(Ka_win_size,Ka_RMS_win,Ka_RMS_win_sigma,'color',c10)
%         hold on
    % plot K1
        plot(K1_win_size,K1_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K1_distalong_transect_m)-min(K1_distalong_transect_m)),K1_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K1_win_size,K1_RMS_win,K1_RMS_win_sigma,'color',c10)
%         hold on
    % plot K2
        plot(K2_win_size,K2_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K2_distalong_transect_m)-min(K2_distalong_transect_m)),K2_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K2_win_size,K2_RMS_win,K2_RMS_win_sigma,'color',c10)
%         hold on
    % plot K4
        plot(K4_win_size,K4_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K4_distalong_transect_m)-min(K4_distalong_transect_m)),K4_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K4_win_size,K4_RMS_win,K4_RMS_win_sigma,'color',c10)
%         hold on
    % plot K5
        plot(K5_win_size,K5_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot((max(K5_distalong_transect_m)-min(K5_distalong_transect_m)),K5_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(K5_win_size,K5_RMS_win,K5_RMS_win_sigma,'color',c10)
%         hold on
    % plot EMC - P, LS
        plot(E1_win_size,E1_RMS_win,'-','color',c10,'linewidth',2)
        hold on
        plot(E1_distalong_transect_m(end)/1000,E1_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E1_win_size,E1_RMS_win,E1_RMS_win_sigma,'color',c10)
%         hold on
    % plot EMC - Borrego
        plot(E2_win_size,E2_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E2_distalong_transect_m(end)/1000,E2_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E2_win_size,E2_RMS_win,E2_RMS_win_sigma,'color',c10)
%         hold on
    % plot EMC - PI
        plot(E3_win_size,E3_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E3_distalong_transect_m(end)/1000,E3_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E3_win_size,E3_RMS_win,E3_RMS_win_sigma,'color',c10) 
%         hold on
    % plot EMC - PS
        plot(E45_win_size,E45_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E45_distalong_transect_m(end)/1000,E45_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E45_win_size,E45_RMS_win,E45_RMS_win_sigma,'color',c10)   
%         hold on
    % plot EMC - Puerta
        plot(E8_win_size,E8_RMS_win,'color',c10,'linewidth',2)
        hold on
        plot(E8_distalong_transect_m(end)/1000,E8_sigma_RMS_dt,'*','color',c10,'MarkerSize',12)
        hold on
%         errorbar(E8_win_size,E8_RMS_win,E8_RMS_win_sigma,'color',c10)
%         hold on
    % plot Balochistan
        plot(B_win_size,B_RMS_win,'-','color',cB,'linewidth',2)
        hold on
        plot(B_distalong_transect_m(end)/1000,B_sigma_RMS_dt,'*','color',cB,'MarkerSize',12)
        hold on
        errorbar(B_win_size,B_RMS_win,B_RMS_win_sigma,'color',cB)
        hold on
    % plot Balochistan - Zinke
        plot(Bz_win_size,Bz_RMS_win,'-','color',cBz,'linewidth',2)
        hold on
        plot(Bz_distalong_transect_m(end)/1000,Bz_sigma_RMS_dt,'*','color',cBz,'MarkerSize',12)
        hold on
        errorbar(Bz_win_size,Bz_RMS_win,Bz_RMS_win_sigma,'color',cBz)
        hold on
    % plot Hector Mine
        plot(HM_win_size,HM_RMS_win,'-','color',cHM,'linewidth',2)
        hold on
        plot(HM_distalong_transect_m(end)/1000,HM_sigma_RMS_dt,'*','color',cHM,'MarkerSize',12)
        hold on
        errorbar(HM_win_size,HM_RMS_win,HM_RMS_win_sigma,'color',cHM)
        hold on
    % plot Landers
        plot(L_win_size,L_RMS_win,'-','color',cL,'linewidth',2)
        hold on
        plot(L_distalong_transect_m(end)/1000,L_sigma_RMS_dt,'*','color',cL,'MarkerSize',12)
        hold on
        errorbar(L_win_size,L_RMS_win,L_RMS_win_sigma,'color',cL)
        hold on
    % plot Palu
        plot(P_win_size,P_RMS_win,'-','color',cP,'linewidth',2)
        hold on
        plot(P_distalong_transect_m(end)/1000,P_sigma_RMS_dt,'*','color',cP,'MarkerSize',12)
        hold on
        errorbar(P_win_size,P_RMS_win,P_RMS_win_sigma,'color',cP)
        hold on
        
         % make legend
            h3(1) = plot(NaN,NaN,'-','color',cB,'linewidth',2);
            h3(2) = plot(NaN,NaN,'-','color',cBz,'linewidth',2);
            h3(3) = plot(NaN,NaN,'-','color',cHM,'linewidth',2);
            h3(4) = plot(NaN,NaN,'-','color',cL,'linewidth',2);
            h3(5) = plot(NaN,NaN,'-','color',cP,'linewidth',2);
            legend(h3,'Balochistan','Balchistan Zinke','Hector Mine','Landers','Palu')
        
            xlim([3*10^-1 3*10^2])
            ylim([7*10^-2 4])
            xlabel('length of window for calculation (km)')
%             ylabel('average value of RMS variance (m)')
            set(gca,'fontsize',16,'fontweight','bold','yscale','log','xscale','log')
            grid on