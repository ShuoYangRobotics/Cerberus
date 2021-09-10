folder = '/home/biorobotics/tightly-coupled-visual-inertial-leg-odometry/output/dataset1_0903/';
% name convention
%  vio - vins-fushion
%  vilo - wob      vilo no bias correction
%  vilo - wb       vilo with bias correction
%  vilo - wb2      vilo with baseline bias correction

% read vio, vilo wob,  vilo wb,
vio_file = strcat(folder,'vio_0903_forward_reidx_with_gt2.csv');
vio_wob_file = strcat(folder,'vio_wlwb_0903_forward_reidx_with_gt2_best.csv');
vio_wb_file = strcat(folder,'vilo_with_leg_bias_0907_best2.csv');

% read table 
vio_Tab = readtable(vio_file);
vio_wob_Tab = readtable(vio_wob_file);
vio_wb_Tab = readtable(vio_wb_file);
data_start_idx = 1;
data_end_idx = 1;

% first read time, find the one with minimum range 
vio_time = (vio_Tab.Var1-vio_Tab.Var1(1))/10^9;
vio_wob_time = (vio_wob_Tab.Var1-vio_wob_Tab.Var1(1))/10^9;
vio_wb_time = (vio_wb_Tab.Var1-vio_wb_Tab.Var1(1))/10^9;

times = {vio_time,vio_wob_time,vio_wb_time};
end_time = [vio_time(end) vio_wob_time(end) vio_wb_time(end)];
[M,I] = min(end_time);

ref_time = times{I};

% read groundtruth position
gt_pos_x = interp1(vio_wb_time,vio_wb_Tab.Var12,ref_time);
gt_pos_y = interp1(vio_wb_time,vio_wb_Tab.Var13,ref_time);
gt_pos_z = interp1(vio_wb_time,vio_wb_Tab.Var14,ref_time);

%calculate ground truth velocity
gt_vel_x = (gt_pos_x(2:end) - gt_pos_x(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
gt_vel_y = (gt_pos_y(2:end) - gt_pos_y(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
gt_vel_z = (gt_pos_z(2:end) - gt_pos_z(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
% smooth the velocity
gt_vel_x = movmean(gt_vel_x,5,1);
gt_vel_y = movmean(gt_vel_y,5,1);
gt_vel_z = movmean(gt_vel_z,5,1);


%% read vio pos
% interpolate time to get all pos
vio_pos_x = interp1(vio_time,vio_Tab.Var2,ref_time);
vio_pos_y = interp1(vio_time,vio_Tab.Var3,ref_time);
vio_pos_z = interp1(vio_time,vio_Tab.Var4,ref_time);
% move data to align with ground truth
vio_pos_x = vio_pos_x + gt_pos_x(1) - vio_pos_x(1);
vio_pos_y = vio_pos_y + gt_pos_y(1) - vio_pos_y(1);
vio_pos_z = vio_pos_z + gt_pos_z(1) - vio_pos_z(1);
% rotate data
angle = -1.9/180*pi;
R = [cos(angle)   -sin(angle)  0;
    sin(angle)    cos(angle)  0;
    0                   0     1];
rotated = R * [vio_pos_x';vio_pos_y';vio_pos_z'];
vio_pos_x = rotated(1,:)';
vio_pos_y = rotated(2,:)';
vio_pos_z = rotated(3,:)';

%calculate vio velocity
vio_vel_x = (vio_pos_x(2:end) - vio_pos_x(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
vio_vel_y = (vio_pos_y(2:end) - vio_pos_y(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
vio_vel_z = (vio_pos_z(2:end) - vio_pos_z(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
% smooth the velocity
vio_vel_x = movmean(vio_vel_x,5,1);
vio_vel_y = movmean(vio_vel_y,5,1);
vio_vel_z = movmean(vio_vel_z,5,1);


%% read vilo wob pos
% interpolate time to get all pos
vio_wob_pos_x = interp1(vio_wob_time,vio_wob_Tab.Var2,ref_time);
vio_wob_pos_y = interp1(vio_wob_time,vio_wob_Tab.Var3,ref_time);
vio_wob_pos_z = interp1(vio_wob_time,vio_wob_Tab.Var4,ref_time);
% move data to align with ground truth
vio_wob_pos_x = vio_wob_pos_x + gt_pos_x(1) - vio_wob_pos_x(1);
vio_wob_pos_y = vio_wob_pos_y + gt_pos_y(1) - vio_wob_pos_y(1);
vio_wob_pos_z = vio_wob_pos_z + gt_pos_z(1) - vio_wob_pos_z(1);
% rotate data
angle = -1.9/180*pi;
R = [cos(angle)   -sin(angle)  0;
    sin(angle)    cos(angle)  0;
    0                   0     1];
rotated = R * [vio_wob_pos_x';vio_wob_pos_y';vio_wob_pos_z'];
vio_wob_pos_x = rotated(1,:)';
vio_wob_pos_y = rotated(2,:)';
vio_wob_pos_z = rotated(3,:)';


%calculate vilo wob velocity
vio_wob_vel_x = (vio_wob_pos_x(2:end) - vio_wob_pos_x(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
vio_wob_vel_y = (vio_wob_pos_y(2:end) - vio_wob_pos_y(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
vio_wob_vel_z = (vio_wob_pos_z(2:end) - vio_wob_pos_z(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
% smooth the velocity
vio_wob_vel_x = movmean(vio_wob_vel_x,5,1);
vio_wob_vel_y = movmean(vio_wob_vel_y,5,1);
vio_wob_vel_z = movmean(vio_wob_vel_z,5,1);

%% read vilo wb pos
% interpolate time to get all pos
vio_wb_pos_x = interp1(vio_wb_time,vio_wb_Tab.Var2,ref_time);
vio_wb_pos_y = interp1(vio_wb_time,vio_wb_Tab.Var3,ref_time);
vio_wb_pos_z = interp1(vio_wb_time,vio_wb_Tab.Var4,ref_time);
% move data to align with ground truth
vio_wb_pos_x = vio_wb_pos_x + gt_pos_x(1) - vio_wb_pos_x(1);
vio_wb_pos_y = vio_wb_pos_y + gt_pos_y(1) - vio_wb_pos_y(1);
vio_wb_pos_z = vio_wb_pos_z + gt_pos_z(1) - vio_wb_pos_z(1);
% rotate data
angle = -1.9/180*pi;
R = [cos(angle)   -sin(angle)  0;
    sin(angle)    cos(angle)  0;
    0                   0     1];
rotated = R * [vio_wb_pos_x';vio_wb_pos_y';vio_wb_pos_z'];
vio_wb_pos_x = rotated(1,:)';
vio_wb_pos_y = rotated(2,:)';
vio_wb_pos_z = rotated(3,:)';

% angle = 3.9/180*pi;
% R = [cos(angle)  0  -sin(angle) ;
%     0  1  0  ;
%     sin(angle)           0        cos(angle)];
% rotated = R * [vio_wb_pos_x;vio_wb_pos_y;vio_wb_pos_z];
% vio_wb_pos_x = rotated(1,:);
% vio_wb_pos_y = rotated(2,:);
% vio_wb_pos_z = rotated(3,:);

%calculate vilo wob velocity
vio_wb_vel_x = (vio_wb_pos_x(2:end) - vio_wb_pos_x(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
vio_wb_vel_y = (vio_wb_pos_y(2:end) - vio_wb_pos_y(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
vio_wb_vel_z = (vio_wb_pos_z(2:end) - vio_wb_pos_z(1:end-1))./(ref_time(2:end)-ref_time(1:end-1));
% smooth the velocity
vio_wb_vel_x = movmean(vio_wb_vel_x,5,1);
vio_wb_vel_y = movmean(vio_wb_vel_y,5,1);
vio_wb_vel_z = movmean(vio_wb_vel_z,5,1);


%% plot all z pos

fontsize = 18
figure(1);clf;hold on;
set(gca,'FontSize',fontsize)
plot(ref_time, gt_pos_z,'LineWidth',2);
plot(ref_time, vio_pos_z,'LineWidth',2);
plot(ref_time, vio_wob_pos_z,'LineWidth',2);
plot(ref_time, vio_wb_pos_z,'LineWidth',2);
legend('Ground Truth', 'VIO', "VILO+NoBiasCorrect", "VILO+BiasCorrect", 'Location','southeast')
xlabel('Time (s)')
ylabel('Z Position (m)')
hold off;
set(gcf, 'Position', [214 1703 550 420])

%% plot all trajectory 
figure(2);clf;
set(gca,'FontSize',fontsize)
plot3(gt_pos_x, gt_pos_y, gt_pos_z,'LineWidth',2);hold on;
plot3(vio_pos_x, vio_pos_y, vio_pos_z,'LineWidth',2);hold on;
plot3(vio_wob_pos_x, vio_wob_pos_y, vio_wob_pos_z,'LineWidth',2);hold on;
plot3(vio_wb_pos_x, vio_wb_pos_y, vio_wb_pos_z,'LineWidth',2);hold off;
legend('Ground Truth', 'VIO', "VILO+NoBiasCorrect", "VILO+BiasCorrect", 'Location','southeast')

xlabel('X Position (m)')
ylabel('Y Position (m)')
zlabel('Z Position (m)')
ax = gca;ax.XLim = [-0.1 3.5];
set(gcf, 'Position', [12 1150 1581 464])


%% plot all x velocity
figure(3);clf;
subplot(3,1,1)
plot(ref_time(1:end-1), gt_vel_x,'LineWidth',2);hold on;
plot(ref_time(1:end-1), vio_vel_x,'LineWidth',2);
plot(ref_time(1:end-1), vio_wob_vel_x,'LineWidth',2);
plot(ref_time(1:end-1), vio_wb_vel_x,'LineWidth',2);
ax = gca;ax.XLim = [20 35];
set(gca,'FontSize',fontsize)
xlabel('Time (s)')
ylabel('X Velocity (m/s)')
legend('Ground Truth', 'VIO', "VILO+NoBiasCorrect", "VILO+BiasCorrect", 'Location','northwest', 'FontSize', 14)
subplot(3,1,2)
plot(ref_time(1:end-1), gt_vel_y,'LineWidth',2);hold on;
plot(ref_time(1:end-1), vio_vel_y,'LineWidth',2);
plot(ref_time(1:end-1), vio_wob_vel_y,'LineWidth',2);
plot(ref_time(1:end-1), vio_wb_vel_y,'LineWidth',2);
set(gca,'FontSize',fontsize)
ax = gca;ax.XLim = [20 35];
xlabel('Time (s)')
ylabel('Y Velocity (m/s)')
legend('Ground Truth', 'VIO', "VILO+NoBiasCorrect", "VILO+BiasCorrect", 'Location','northwest', 'FontSize', 14)
subplot(3,1,3)
plot(ref_time(1:end-1), gt_vel_z,'LineWidth',2);hold on;
plot(ref_time(1:end-1), vio_vel_z,'LineWidth',2);
plot(ref_time(1:end-1), vio_wob_vel_z,'LineWidth',2);
plot(ref_time(1:end-1), vio_wb_vel_z,'LineWidth',2);
set(gca,'FontSize',fontsize)
ax = gca;ax.XLim = [20 35];
xlabel('Time (s)')
ylabel('Z Velocity (m/s)')
legend('Ground Truth', 'VIO', "VILO+NoBiasCorrect", "VILO+BiasCorrect", 'Location','northwest', 'FontSize', 14)
set(gcf, 'Position', [20 27 1567 1038])