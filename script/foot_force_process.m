bagselect = rosbag('best_x.bag');
bSel = select(bagselect,'Topic','/hardware_a1/joint_foot');
msgStructs = readMessages(bSel,'DataFormat','struct');

num_data = size(msgStructs,1);
num_leg = 4;

start_idx = 3000;
end_idx = 5110;

foot_force = zeros(num_data,4);
for i=1:num_data
    foot_force(i,1) = msgStructs{i}.Effort(13);
    foot_force(i,2) = msgStructs{i}.Effort(14);
    foot_force(i,3) = msgStructs{i}.Effort(15);
    foot_force(i,4) = msgStructs{i}.Effort(16);
end


foot_force_extrema = zeros(4,2);
foot_force_min = zeros(num_data,4);
foot_force_max = zeros(num_data,4);

% calculate contact flag
contact_flag = zeros(num_data,4);
threadholds = zeros(num_data,4);
for i=1:num_data
    for j=1:num_leg
        if (foot_force(i,j) < foot_force_extrema(j,1))
            foot_force_extrema(j,1) = foot_force_extrema(j,1) * 0.9;
            
            foot_force_extrema(j,1) = foot_force_extrema(j,1) + 0.1*foot_force(i,j);
        end
        if (foot_force(i,j) > foot_force_extrema(j,2))
            foot_force_extrema(j,2) = foot_force_extrema(j,2) * 0.9;
            
            foot_force_extrema(j,2) = foot_force_extrema(j,2) + 0.1*foot_force(i,j);
        end
        foot_force_extrema(j,1) = foot_force_extrema(j,1) * 0.9991;
        foot_force_extrema(j,2) = foot_force_extrema(j,2) * 0.997;
        
        % record extreme values
        foot_force_min(i,j) = foot_force_extrema(j,1);
        foot_force_max(i,j) = foot_force_extrema(j,2);

        threadhold = 0.4*(foot_force_extrema(j,2)-foot_force_extrema(j,1))+foot_force_extrema(j,1);
        if (foot_force(i,j) > threadhold)
            contact_flag(i,j) = 1;
        end
        threadholds(i,j) = threadhold;
    end
end



figure(1)
for j=1:1
%     subplot(4,2,(j-1)*2+1)
%     subplot(4,1,(j-1)*1+1)
    plot(foot_force(start_idx:end_idx,j),'b', "LineWidth",4);hold on;
    %plot(foot_force_min(start_idx:end_idx,j),'r');hold on;
    %plot(foot_force_max(start_idx:end_idx,j),'r');hold on;
%     subplot(4,2,(j-1)*2+2)
    plot(threadholds(start_idx:end_idx,j),'y', "LineWidth",4);hold on;
    plot(150*contact_flag(start_idx:end_idx,j),'r', "LineWidth",4);hold off;
    legend("Sensor Reading","Threshold","Contact flag (Scaled)")
    set(gca,'FontSize',18)
end