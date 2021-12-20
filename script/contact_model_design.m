% the min and max of the noise
n_min = 0.2;
n_max = 2000;
% weight
w1 = 0.3;
w3 = 0.6;
w2 = 1-w1-w3;
% term 1 parameters
steep = 1;
threshold = 50;

% term 2 parameters
dt = 1/500;
bound = 1000/dt %in one tick the force changes 500 
rescale = 0.3;

% term 3 parameters
distance_scale = 1000;

% final adjust
final_ratio = 0.1;

% sensor input 
% previous and current force
c_0 = 250;
c_1 = 249;
% the velocity of the base
base_v = [0.4;0.0;0.1];

% the calculated leg velocity
lo_v = [1.7;0.1;0.1]; 

% term 1 foot force thresholding 
force_mag = 0.5*(c_0+c_1);
n1 = n_max*(1-1/(1+exp(-steep*(force_mag-threshold))))+n_min

% term 2 variation 
diff_force = (c_1 - c_0)/dt;
diff_force = diff_force/bound;
if (diff_force > 0) 
    diff_force = 0.5*diff_force;
end
n2 = abs(diff_force)*n_max

% test ( assume we load a bag using foot_force_process.m
diff_list=(foot_force(2:end,1)-foot_force(1:end-1,1))/dt;
diff_list = diff_list/bound;
diff_list(diff_list>0) = rescale*diff_list(diff_list>0);
diff_list = abs(diff_list);
diff_list = max(0,sqrt(diff_list)*n_max-30);
plot(diff_list)

% term 3 
diff_v = lo_v - base_v;
n3 = distance_scale*max(0,sqrt((diff_v).^2) - [0.2;0.2;0.1])+n_min;

n = final_ratio*(w1*n1*ones(3,1)+w2*n2*ones(3,1)+w3*n3)