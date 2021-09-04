folder = 'I:\research4year\tightly-coupled-visual-inertial-leg-odometry\output\';
% read vio
file = strcat(folder,'vio_with_leg_without_bias_estimation015.csv');
T = readtable(file);
twwo015 = (T.Var1-T.Var1(1))/10^9;
vio_wwo015_pos_x = T.Var2;
vio_wwo015_pos_y = T.Var3;
vio_wwo015_pos_z = T.Var4;

% plot(twwo015, vio_wwo015_pos_z,'LineWidth',2);

% read vio with_leg_without_bias_correction
file = strcat(folder,'vio_with_leg_without_bias_estimation020.csv');
T = readtable(file);
twwo020 = (T.Var1-T.Var1(1))/10^9;
vio_wwo020_pos_x = T.Var2;
vio_wwo020_pos_y = T.Var3;
vio_wwo020_pos_z = T.Var4;

plot(twwo020, vio_wwo020_pos_z,'LineWidth',2);hold on;

% read vio with_leg_without_bias_correction
file = strcat(folder,'vio_with_leg_without_bias_estimation025.csv');
T = readtable(file);
twwo025 = (T.Var1-T.Var1(1))/10^9;
vio_wwo025_pos_x = T.Var2;
vio_wwo025_pos_y = T.Var3;
vio_wwo025_pos_z = T.Var4;

plot(twwo025, vio_wwo025_pos_z,'LineWidth',2)

% legend("VIO+Leg/No bias correction, Lm=0.15m", "VIO+Leg/No bias correction, Lm=0.20m", "VIO+Leg/No bias correction, Lm=0.25m", 'Location','southeast')
legend( "VIO+Leg/No bias correction, Lm=0.20m", "VIO+Leg/No bias correction, Lm=0.25m", 'Location','southeast')
% read vio
file = strcat(folder,'vio_with_leg_with_bias_estimation015.csv');
T = readtable(file);
tww015 = (T.Var1-T.Var1(1))/10^9;
vio_ww015_pos_x = T.Var2;
vio_ww015_pos_y = T.Var3;
vio_ww015_pos_z = T.Var4;
figure
% plot(tww015, vio_ww015_pos_z,'LineWidth',2);hold on;

% read vio with_leg_without_bias_correction
file = strcat(folder,'vio_with_leg_with_bias_estimation020.csv');
T = readtable(file);
tww020 = (T.Var1-T.Var1(1))/10^9;
vio_ww020_pos_x = T.Var2;
vio_ww020_pos_y = T.Var3;
vio_ww020_pos_z = T.Var4;

plot(tww020, vio_ww020_pos_z,'LineWidth',2);;hold on;

% read vio with_leg_without_bias_correction
file = strcat(folder,'vio_with_leg_with_bias_estimation025.csv');
T = readtable(file);
tww025 = (T.Var1-T.Var1(1))/10^9;
vio_ww025_pos_x = T.Var2;
vio_ww025_pos_y = T.Var3;
vio_ww025_pos_z = T.Var4;

plot(tww025, vio_ww025_pos_z,'LineWidth',2)

% plot(t3, 0.21*ones(size(t3)),'--','LineWidth',2)
% legend("VIO+Leg/With bias correction, Lm=0.15m", "VIO+Leg/With bias correction, Lm=0.20m", "VIO+Leg/With bias correction, Lm=0.25m", 'Location','southeast')

legend("VIO+Leg/With bias correction, Lm=0.20m", "VIO+Leg/With bias correction, Lm=0.25m", 'Location','southeast')
% legend('VIO', "VIO+Leg/No bias correction", "VIO+Leg/With bias correction", "Stance height", 'Location','southeast')

% trajectory plot
figure 

plot3(vio_wwo015_pos_x, vio_wwo015_pos_y, vio_wwo015_pos_z,'LineWidth',2);hold on;
plot3(vio_wwo020_pos_x, vio_wwo020_pos_y, vio_wwo020_pos_z,'LineWidth',2);hold on;
plot3(vio_wwo025_pos_x, vio_wwo025_pos_y, vio_wwo025_pos_z,'LineWidth',2);hold on;
plot3(vio_ww015_pos_x, vio_ww015_pos_y, vio_ww015_pos_z,'LineWidth',2);hold on;
plot3(vio_ww020_pos_x, vio_ww020_pos_y, vio_ww020_pos_z,'LineWidth',2);hold on;
plot3(vio_ww025_pos_x, vio_ww025_pos_y, vio_ww025_pos_z,'LineWidth',2);hold on;
axis equal

legend("VIO+Leg/No bias correction, Lm=0.15m", "VIO+Leg/No bias correction, Lm=0.20m", "VIO+Leg/No bias correction, Lm=0.25m",...
    "VIO+Leg/With bias correction, Lm=0.15m", "VIO+Leg/With bias correction, Lm=0.20m", "VIO+Leg/With bias correction, Lm=0.25m", 'Location','southeast')
