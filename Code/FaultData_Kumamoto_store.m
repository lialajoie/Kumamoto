% data plotting limits
% FullData = [617550 662550];
% FullData = [3570050 3607450];
% Fault sections: 1:77

buffer = 1000; % amount of space in x and y direction away from fault

%% FAULT ALL
FaultData(10).FaultName = 'Kumamoto';
FaultData(10).SegmentName = 'ALL, 45 azimuth'; 
FaultData(10).ZDispRange = [-2 0.5];
FaultData(10).XDispRange = [-2 2];
FaultData(10).Xlim = [-23856 -7376];
FaultData(10).Ylim = [-31368 -17376];
FaultData(10).Xfault = [-8754 -21000];
FaultData(10).Yfault = [-17500 -30000];

%% FAULT ALL - GENERALIZED WITH THREE SEGMENTS
FaultData(11).FaultName = 'Kumamoto - ALL';
FaultData(11).SegmentName = 'Generalized, 3 segments'; 
FaultData(11).ZDispRange = [-2 0.5];
FaultData(11).XDispRange = [-2 2];
FaultData(11).map_Xlim = [-23856 -7376];
FaultData(11).map_Ylim = [-31368 -17376];
FaultData(11).Xfault = [-8544 -11760 -18600 -19472]; 
FaultData(11).Yfault = [-19314 -20496 -27000 -30312];

%% FAULT 1
X1 = [-20127 -18652];
Y1 = [-31932 -27818];
FaultData(1).FaultName = 'Southern Zone';
FaultData(1).SegmentName = 'Southern segment'; 
FaultData(1).ZDispRange = [-2 0.5];
FaultData(1).XDispRange = [-2 2];
FaultData(1).map_Xlim = [X1(1)-buffer X1(2)+buffer];
FaultData(1).map_Ylim = [Y1(1)-buffer Y1(2)+buffer];
FaultData(1).Xfault = [-19428 -18837];
FaultData(1).Yfault = [-30121 -27818];

%% FAULT 2
X2 = [-17614 -15551];
Y2 = [-26302 -24681];
FaultData(2).FaultName = 'Central Zone';
FaultData(2).SegmentName = 'Southern segment'; 
FaultData(2).ZDispRange = [-2 0.5];
FaultData(2).XDispRange = [-2 2];
FaultData(2).map_Xlim = [X2(1)-buffer X2(2)+buffer];
FaultData(2).map_Ylim = [Y2(1)-buffer Y2(2)+buffer];
FaultData(2).Xfault = [-17614 -15551];
FaultData(2).Yfault = [-26302 -24681];

%% FAULT 3
X3 = [-16215 -13525];
Y3 = [-24456 -22336];
FaultData(3).FaultName = 'Central Zone';
FaultData(3).SegmentName = 'Central segment'; 
FaultData(3).ZDispRange = [-2 0.5];
FaultData(3).XDispRange = [-2 2];
FaultData(3).map_Xlim = [X3(1)-buffer X3(2)+buffer];
FaultData(3).map_Ylim = [Y3(1)-buffer Y3(2)+buffer];
FaultData(3).Xfault = [-16215 -13525];
FaultData(3).Yfault = [-24456 -22336];

%% FAULT 4
X4 = [-14164 -13744];
Y4 = [-22431 -22213];
FaultData(4).FaultName = 'Central Zone';
FaultData(4).SegmentName = 'Oblique connection'; 
FaultData(4).ZDispRange = [-2 0.5];
FaultData(4).XDispRange = [-2 2];
FaultData(4).map_Xlim = [X4(1)-buffer X4(2)+buffer];
FaultData(4).map_Ylim = [Y4(1)-buffer Y4(2)+buffer];
FaultData(4).Xfault = [-13747 -14232];
FaultData(4).Yfault = [-22430 -22168];

%% FAULT 5
X5 = [-16165 -12515];
Y5 = [-23084 -21120];
FaultData(5).FaultName = 'Central Zone';
FaultData(5).SegmentName = 'Northern segment'; 
FaultData(5).ZDispRange = [-2 0.5];
FaultData(5).XDispRange = [-2 2];
FaultData(5).map_Xlim = [X5(1)-buffer X5(2)+buffer];
FaultData(5).map_Ylim = [Y5(1)-buffer Y5(2)+buffer];
FaultData(5).Xfault = [-16164 -12515];
FaultData(5).Yfault = [-23084 -21120];

%% FAULT 6
X6 = [-9710 1183];
Y6 = [-19868 -11298];
FaultData(6).FaultName = 'Northern Zone';
FaultData(6).SegmentName = 'Northern segment'; 
FaultData(6).ZDispRange = [-2 0.5];
FaultData(6).XDispRange = [-2 2];
FaultData(6).map_Xlim = [X6(1)-buffer X6(2)+buffer];
FaultData(6).map_Ylim = [Y6(1)-buffer Y6(2)+buffer];
FaultData(6).Xfault = [];
FaultData(6).Yfault = [];

%% FAULT 7
X7 = [-7684 -4272];
Y7 = [-19973 -17852];
FaultData(7).FaultName = 'Northern Zone';
FaultData(7).SegmentName = 'Southern segment'; 
FaultData(7).ZDispRange = [-2 0.5];
FaultData(7).XDispRange = [-2 2];
FaultData(7).map_Xlim = [X7(1)-buffer X7(2)+buffer];
FaultData(7).map_Ylim = [Y7(1)-buffer Y7(2)+buffer];
FaultData(7).Xfault = [];
FaultData(7).Yfault = [];







