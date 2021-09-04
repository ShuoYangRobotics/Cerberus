folder = 'I:\research4year\tightly-coupled-visual-inertial-leg-odometry\output\';
% read vio
file = strcat(folder,'vio.csv');
T = readtable(file);
t1 = (T.Var1-T.Var1(1))/10^9;
vio_pos_x = T.Var2;
vio_pos_y = T.Var3;
vio_pos_z = T.Var4;

plot(t1, vio_pos_z,'LineWidth',2);hold on;

% read vio with_leg_without_bias_correction
file = strcat(folder,'vio_with_leg_without_bias_estimation.csv');
T = readtable(file);
t2 = (T.Var1-T.Var1(1))/10^9;
vio_wwo_pos_x = T.Var2;
vio_wwo_pos_y = T.Var3;
vio_wwo_pos_z = T.Var4;

plot(t2, vio_wwo_pos_z,'LineWidth',2)

% read vio with_leg_with_bias_correction
file = strcat(folder,'vio_with_leg_with_bias_estimation.csv');
T = readtable(file);
t3 = (T.Var1-T.Var1(1))/10^9;
vio_ww_pos_x = T.Var2;
vio_ww_pos_y = T.Var3;
vio_ww_pos_z = T.Var4;

% height drift comparison
plot(t3, vio_ww_pos_z,'LineWidth',2)

% plot(t3, 0.21*ones(size(t3)),'--','LineWidth',2)
legend('VIO', "VIO+Leg/No bias correction", "VIO+Leg/With bias correction", "Stance height", 'Location','southeast')
% legend('VIO', "VIO+Leg/No bias correction", "VIO+Leg/With bias correction", "Stance height", 'Location','southeast')

% trajectory plot
figure 

plot3(vio_pos_x, vio_pos_y, vio_pos_z,'LineWidth',2);hold on;
plot3(vio_wwo_pos_x, vio_wwo_pos_y, vio_wwo_pos_z,'LineWidth',2);hold on;
plot3(vio_ww_pos_x, vio_ww_pos_y, vio_ww_pos_z,'LineWidth',2);hold on;
axis equal

legend('VIO', "VIO+Leg/No bias correction", "VIO+Leg/With bias correction",'Location','southeast')
