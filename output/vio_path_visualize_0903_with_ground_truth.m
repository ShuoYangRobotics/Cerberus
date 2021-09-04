folder = 'I:\research4year\tightly-coupled-visual-inertial-leg-odometry\output\';
% read vio
% file = strcat(folder,'vio_0903_forward_reidx_with_gt.csv');
data_start_idx = 100;
file = strcat(folder,'vio_0903_forward_backward_reidx_with_gt.csv');
T = readtable(file);
data_end_idx = size(T.Var1,1)-50;
t1 = (T.Var1(data_start_idx:data_end_idx)-T.Var1(data_start_idx))/10^9;
vio_pos_x = T.Var2(data_start_idx:data_end_idx);
vio_pos_y = T.Var3(data_start_idx:data_end_idx);
vio_pos_z = T.Var4(data_start_idx:data_end_idx);

% move gt to align with vio
offset_x =  T.Var2(data_start_idx) - T.Var12(data_start_idx);
offset_y =  T.Var3(data_start_idx) - T.Var13(data_start_idx);
offset_z =  T.Var4(data_start_idx) - T.Var14(data_start_idx);

vio_gt_pos_x = T.Var12(data_start_idx:data_end_idx) + offset_x;
vio_gt_pos_y = T.Var13(data_start_idx:data_end_idx) + offset_y;
vio_gt_pos_z = T.Var14(data_start_idx:data_end_idx) + offset_z;

plot(t1, vio_gt_pos_z,'LineWidth',2);hold on;
plot(t1, vio_pos_z,'LineWidth',2);hold on;

% read vio with leg with _bias_correction
% file = strcat(folder,'vio_wlwb_0903_forward_reidx_with_gt.csv');
file = strcat(folder,'vio_wlwb_0903_forward_backward_reidx_with_gt.csv');
T = readtable(file);
data_end_idx = size(T.Var1,1)-100;
t_wlwb = (T.Var1(data_start_idx:data_end_idx)-T.Var1(data_start_idx))/10^9;
vio_wlwb_pos_x = T.Var2(data_start_idx:data_end_idx);
vio_wlwb_pos_y = T.Var3(data_start_idx:data_end_idx);
vio_wlwb_pos_z = T.Var4(data_start_idx:data_end_idx);

% move gt to align with vio
offset_x =  T.Var2(data_start_idx) - T.Var12(data_start_idx)
offset_y =  T.Var3(data_start_idx) - T.Var13(data_start_idx)
offset_z =  T.Var4(data_start_idx) - T.Var14(data_start_idx)

vio_wlwb_gt_pos_x = T.Var12(data_start_idx:data_end_idx) + offset_x;
vio_wlwb_gt_pos_y = T.Var13(data_start_idx:data_end_idx) + offset_y;
vio_wlwb_gt_pos_z = T.Var14(data_start_idx:data_end_idx) + offset_z;

plot(t_wlwb, vio_wlwb_pos_z,'LineWidth',2);hold on;
legend('Ground Truth', 'VIO', "VILO+Proposed",'Location','southeast')
% plot(t_wlwb, vio_wlwb_gt_pos_z,'LineWidth',2);hold on;

% % read vio with_leg_without_bias_correction
% file = strcat(folder,'vio_with_leg_without_bias_estimation.csv');
% T = readtable(file);
% t2 = (T.Var1-T.Var1(1))/10^9;
% vio_wwo_pos_x = T.Var2;
% vio_wwo_pos_y = T.Var3;
% vio_wwo_pos_z = T.Var4;
% 
% plot(t2, vio_wwo_pos_z,'LineWidth',2)
% 
% % read vio with_leg_with_bias_correction
% file = strcat(folder,'vio_with_leg_with_bias_estimation.csv');
% T = readtable(file);
% t3 = (T.Var1-T.Var1(1))/10^9;
% vio_ww_pos_x = T.Var2;
% vio_ww_pos_y = T.Var3;
% vio_ww_pos_z = T.Var4;
% 
% % height drift comparison
% plot(t3, vio_ww_pos_z,'LineWidth',2)
% 
% % plot(t3, 0.21*ones(size(t3)),'--','LineWidth',2)
% legend('VIO', "VIO+Leg/No bias correction", "VIO+Leg/With bias correction", "Stance height", 'Location','southeast')
% % legend('VIO', "VIO+Leg/No bias correction", "VIO+Leg/With bias correction", "Stance height", 'Location','southeast')
% 
% trajectory plot
figure 

plot3(vio_gt_pos_x, vio_gt_pos_y, vio_gt_pos_z,'LineWidth',2);hold on;
plot3(vio_pos_x, vio_pos_y, vio_pos_z,'LineWidth',2);hold on;
plot3(vio_wlwb_pos_x, vio_wlwb_pos_y, vio_wlwb_pos_z,'LineWidth',2);hold on;
axis equal

legend('Ground Truth', 'VIO', "VILO+Proposed",'Location','southeast')
