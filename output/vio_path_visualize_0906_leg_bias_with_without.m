folder = '/home/shuoy/nv_work/unitree_ws/src/tightly-coupled-visual-inertial-leg-odometry/output/';
% read vio
file = strcat(folder,'vio_0903_forward_reidx_with_gt2.csv');
data_start_idx = 1;
% file = strcat(folder,'vio_0903_forward_backward_reidx_with_gt.csv');
T = readtable(file);
data_end_idx = size(T.Var1,1)-1;
t1 = (T.Var1(data_start_idx:data_end_idx)-T.Var1(data_start_idx))/10^9;
vio_pos_x = T.Var2(data_start_idx:data_end_idx);
vio_pos_y = T.Var3(data_start_idx:data_end_idx);
vio_pos_z = T.Var4(data_start_idx:data_end_idx);

figure(1)
angle = -1.9/180*pi;
R = [cos(angle)   -sin(angle)  0;
    sin(angle)    cos(angle)  0;
    0                   0     1];
rotated = R * [vio_pos_x';vio_pos_y';vio_pos_z'];
vio_pos_x = rotated(1,:);
vio_pos_y = rotated(2,:);
vio_pos_z = rotated(3,:);


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
file = strcat(folder,'vilo_with_leg_bias_0907_best2.csv');
% file = strcat(folder,'vio_wlwb_0903_forward_backward_reidx_with_gt2.csv');
T = readtable(file);
data_end_idx = size(T.Var1,1)-1;
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

R = [cos(angle)   -sin(angle)  0;
    sin(angle)    cos(angle)  0;
    0                   0     1];
rotated = R * [vio_wlwb_pos_x';vio_wlwb_pos_y';vio_wlwb_pos_z'];
vio_wlwb_pos_x = rotated(1,:);
vio_wlwb_pos_y = rotated(2,:);
vio_wlwb_pos_z = rotated(3,:);

plot(t_wlwb, vio_wlwb_pos_z,'LineWidth',2);hold on;

% read vio with leg with _bias_correction
file = strcat(folder,'vio.csv');
% file = strcat(folder,'vio_wlwb_0903_forward_backward_reidx_with_gt2.csv');
T = readtable(file);
data_end_idx = size(T.Var1,1)-1;
t_wlwob = (T.Var1(data_start_idx:data_end_idx)-T.Var1(data_start_idx))/10^9;
vio_wlwob_pos_x = T.Var2(data_start_idx:data_end_idx);
vio_wlwob_pos_y = T.Var3(data_start_idx:data_end_idx);
vio_wlwob_pos_z = T.Var4(data_start_idx:data_end_idx);

% move gt to align with vio
offset_x =  T.Var2(data_start_idx) - T.Var12(data_start_idx)
offset_y =  T.Var3(data_start_idx) - T.Var13(data_start_idx)
offset_z =  T.Var4(data_start_idx) - T.Var14(data_start_idx)

vio_wlwob_gt_pos_x = T.Var12(data_start_idx:data_end_idx) + offset_x;
vio_wlwob_gt_pos_y = T.Var13(data_start_idx:data_end_idx) + offset_y;
vio_wlwob_gt_pos_z = T.Var14(data_start_idx:data_end_idx) + offset_z;

R = [cos(angle)   -sin(angle)  0;
    sin(angle)    cos(angle)  0;
    0                   0     1];
rotated = R * [vio_wlwob_pos_x';vio_wlwob_pos_y';vio_wlwob_pos_z'];
vio_wlwob_pos_x = rotated(1,:);
vio_wlwob_pos_y = rotated(2,:);
vio_wlwob_pos_z = rotated(3,:);

plot(t_wlwob, vio_wlwob_pos_z,'LineWidth',2);hold off;

legend('Ground Truth', 'VIO', "VILO+Proposed", "VILO+NoBiasCorrect",'Location','southeast')

figure(2)

plot3(vio_gt_pos_x, vio_gt_pos_y, vio_gt_pos_z,'LineWidth',2);hold on;
plot3(vio_pos_x, vio_pos_y, vio_pos_z,'LineWidth',2);hold on;
plot3(vio_wlwb_pos_x, vio_wlwb_pos_y, vio_wlwb_pos_z,'LineWidth',2);hold on;
plot3(vio_wlwob_pos_x, vio_wlwob_pos_y, vio_wlwob_pos_z,'LineWidth',2);hold off;
axis equal

legend('Ground Truth', 'VIO', "VILO+Proposed", "VILO+NoBiasCorrect",'Location','southeast')


% pose diff
% need to interpolate time 
gt_pos = [vio_gt_pos_x(end);vio_gt_pos_y(end);vio_gt_pos_z(end)];
vio_ww_pos = [vio_pos_x(end);vio_pos_y(end);vio_pos_z(end)];
vio_wlwb_pos = [vio_wlwb_pos_x(end);vio_wlwb_pos_y(end);vio_wlwb_pos_z(end)];
norm(vio_ww_pos - gt_pos)
norm(vio_wlwb_pos - gt_pos)
