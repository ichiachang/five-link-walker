% Run Five-link robot walking on DRS simulation
% with one-step ahead controller
% ode45 simulation

clear;
close all;
clc;

% add path ...
addpath("gen_byFROST\")
addpath("sub_func\")

%% DRS motion setting
DRSsetting.ampX = 0.15 ; % m
DRSsetting.ampY = 0;
DRSsetting.ampZ = 0;
DRSsetting.ampRX = 0;
DRSsetting.ampRY = 0;
DRSsetting.ampRZ = 0;

DRSsetting.freqX = 1; % rad/sec
DRSsetting.freqY = 1;
DRSsetting.freqZ = 1;
DRSsetting.freqRX = 1;
DRSsetting.freqRY = 1;
DRSsetting.freqRZ = 1;

% Define the DRS motion as a function handle
DRSmotion_h = @(t)DRSmotion(t,DRSsetting);
DRSmotion2_h = @(t,r,c)DRSmotion2(t,DRSsetting,r,c);
% check statespace function
% x_dot = FiveLink_constrained_dynamics(t0,xR_ini,uR_ini,DRSmotion_h)

%% Simulation setting
T_total_sim = 1; % sec, Duration of simulation
freq_ctl_l = 1000; % Hz, low level control (I/O linearization) freq.
freq_ctl_h = 200;  % Hz, high level control (ALIP) freq.

theta_b_star = 10; % deg
q1_ini = 130; % deg
q2_ini = 80; % deg

%% Initial condition setting

% qR: robot state (Full order) w.r.t. {W} frame
% qR = [x, z, theta, q1_r, q2_r, q1_l, q2_l]' [7x1]
qR_ini = deg2rad([0, 0, theta_b_star, q1_ini, q2_ini, q1_ini, q2_ini]');  % x, z will be solved based on contact constraint
right_toe_pos = pos_toe_r_func(qR_ini);         % such that the right(stance) toe locates at (0,0) w.r.t. {W}
qR_ini(1:2,1) = -[right_toe_pos(1); right_toe_pos(3)];
% % show the robot configuration
figure;
plot_FiveLinkWalker(qR_ini); plot([-1 0 1],[0 0 0],'kx-')

% set up the qR_dot_ini
qR_dot_ini = zeros(7,1);
DRS = DRSmotion(0,DRSsetting);
qR_dot_ini(1) = DRS(1,2); % make the CoM have same speed as DRS

% xR = [ qR; qR_dot ] [14x1]
xR_ini = [qR_ini; qR_dot_ini];
% uR: robot joint input 
% uR = [tau_1r, tau_2r, tau_1l, tau_2l]' [4x1]
uR_ini = zeros(4,1);

%% Robot/Controller parameter setting
Robot.m = 32; % mass, kg, from URDF file
CoM_ini = pos_CoM_func(qR_ini);
Robot.H = CoM_ini(3); % height of robot CoM, height of ALIP model
Robot.g = 9.81; % m/sec^2

controller.T = 0.35; % stepping duration, sec
controller.v_d = 0.2; % desired velocity, m/sec
controller.L_s_d = Robot.m * Robot.H * controller.v_d; % desired angular velocity

controller.Bezier.z = Robot.H * ones(7,1);
controller.Bezier.theta_b = deg2rad(theta_b_star) * ones(7,1);
controller.Bezier.x_sw_base = [0 0.05 0.15 0.50 0.80 0.95 1]'; %[-1 -0.9 -0.7 0 0.7 0.9 1]' + 1;
controller.Bezier.x_sw = controller.Bezier.x_sw_base * controller.T * controller.v_d; 
controller.Bezier.z_sw = [0 0.06 0.07 0.08 0.06 0.03 0]';
controller.BezierAlpha = [controller.Bezier.z';
                          controller.Bezier.theta_b';
                          controller.Bezier.x_sw';
                          controller.Bezier.z_sw']'; % [7x4]

controller.IO_kp = 2500; % controller gain for the input/output linearization controller
controller.IO_kd = 100;  


%% Simulation loop

%%% Variables initialization
xR = xR_ini;
uR = uR_ini;
t_sim_f = 0;
t_last_td = t_sim_f; % previous touchdown timing;
t_next_td = t_last_td + controller.T; % up-coming touchdown timing
x_sw_td = 0;   % when touchdown, the swing foot position w.r.t. world frame
x_st_td = 0;   % when touchdown, the stance foot position w.r.t. world frame
x_drs_td = DRSmotion2_h(0,1,1);  % when touchdown, the DRS position w.r.t. world frame

t_list = [];
x_list = [];
u_list = [];
step_list = [];

t2_list2 = []; % time stamp for each high-level control loop
x_list2 = [];
u_list2 = [];
qA_list2 = []; % ALIP states for each high-level control loop
step_length_list2 = []; % desired step length for each high-level control loop
hd_list2 = []; % desired task space list aligned with t2
hc_list2 = []; % current task space list aligned with t2
hd_cur_list2 = []; % 
step_list2 = [];

step_index = 1;
%%% END


%%%%%%%%%%%%%%%%%%
%%% WHILE LOOP %%%
%%%%%%%%%%%%%%%%%%
isend = 0; % flag for while loop
while isend == 0

    % initialize the simulation time span
    t_sim_i = t_sim_f;
    t_sim_f = t_sim_i + 1/freq_ctl_l ;
    t_sim_span = [t_sim_i, t_sim_f];
    disp(['simulating from ',num2str(t_sim_i),' to ',num2str(t_sim_f)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% HIGH LEVEL: ALIP controller %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 2: Desired task space trajectory generator
    % desired step length by one-step ahead controller
    qA = Robot2ALIP_state(xR);
    t_cur = t_sim_i; %
    step_length = one_step_ahead(qA,t_cur,t_next_td,DRSmotion2_h,Robot,controller);
    % modified the desired task space reference trajectory ( with no
    % consideration of DRS motion )
    [hd,alpha_d] = desired_task_space_traj(step_length,x_sw_td,x_st_td,controller);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LOW LEVEL: I/O linearizatoin %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 3: Input/output linearization controller
    hc = Robot2task_state(xR); % or hc = hc_func(xR);
    
    s_cur = mod(t_cur - t_last_td , controller.T) / controller.T; % calculate the current phase variable, s_cur
    [hd_cur2, dhd_ds_cur2, ddhd_dss_cur2] = Bezier_6th(alpha_d,s_cur);
    hd_delta = [0; 0; DRSmotion2_h(t_cur,1,1) - x_drs_td; 0]; % shift the hd according to the DRS motion
    hd_cur = hd(s_cur) + hd_delta; % [4x1] current desired task space position w.r.t. world frame
    dhd_ds_delta = [0; 0; DRSmotion2_h(t_cur,1,2); 0];
    dhd_ds_cur = dhd_ds_cur2 + dhd_ds_delta; % [4x1] current first time derivative (w.r.t. s) of desired task space position w.r.t. world frame 
    ddhd_dss_delta = [0; 0; DRSmotion2_h(t_cur,1,3); 0];
    ddhd_dss_cur = ddhd_dss_cur2 + ddhd_dss_delta; % [4x1] current sencod time derivative (w.r.t. s) of desired task space position w.r.t. world frame 
    u_IO = IO_linearization(hc,hd_cur,dhd_ds_cur,ddhd_dss_cur,xR,t_cur,DRSsetting,controller);
    % u_IO = 0*u_IO;
    
    %%% 4: Solve ODE with consistant control input within the same control period
    % define switch event
    options=odeset('Events',@(t,x)SE_touchdown(t,x,DRSmotion_h));
    % ode45 to solve the statespace function:
    [t_ode,x_ode,t_td,x_td,i_td] = ode45(@(t,x)FiveLink_constrained_dynamics(t,x,u_IO,DRSmotion_h),...
                                        t_sim_span,xR,options); % t_td, x_td, i_td: the info when touchdown
    
    %%% check the state
    
    % figure;
    % plot_FiveLinkWalker(x_ode(1,1:7));
    % plot_FiveLinkWalker(x_ode(end,1:7));
    
    %%%%%%%%%%%%%%%%%%%%%
    
    %%% 5: record the state
    t_list = [t_list, t_ode(1:end-1,1).'];
    x_list = [x_list, x_ode(1:end-1,:).'];
    u_list = [u_list, u_IO*ones(1,length(t_ode)-1)];
    step_list = [step_list, ones(1,length(t_ode)-1)*step_index];
    
    t2_list2 = [t2_list2, t_sim_i];
    x_list2 = [x_list2, x_ode(1,:).'];
    u_list2 = [u_list2, u_IO];
    qA_list2 = [qA_list2, qA];
    step_length_list2 = [step_length_list2, step_length];
    hc_list2 = [hc_list2, hc];
    hd_list2 = [hd_list2, hd(s_cur)];
    hd_cur_list2 = [hd_cur_list2, hd_cur];
    step_list2 = [step_list2, step_index];
    %%% 5_end: finish recording
    
    %%% 6: Update state variables for each loop
    t_sim_i = t_ode(end,1);
    xR = x_ode(end,:)';
    
    %%% 6+: check whether the touchdown happen, if the touchdown happen before
    %%% the schedule timing, we, here, just finish this step and begin next
    %%% step. 
    if isempty(t_td) == 0
        step_index = step_index + 1;

        if t_ode(end) > T_total_sim
            isend = 1;
        end

        % update important timing
        t_sim_f = t_ode(end); % next control loop starts at the TD event
        t_last_td = t_sim_f;  % update the last TD timing (i.e. now)
        t_next_td = t_last_td + controller.T; % next touchdown timing is the current one plus desired stepping duration
    
        % update swing foot initial position for next step
        x_st_pre_td = pos_toe_r_func(xR(1:7)); % obtain the pre-TD stance(right) toe position
        x_st_pre_td = x_st_pre_td(1);
        x_sw_td = x_st_pre_td; % <<< update the initial swing foot position for next step's Bezier curve
        % the swing(left) leg for next step is the stance(right) leg for this
        % step.
    
        % update stance foot initial position for next step
        x_sw_pre_td = pos_toe_l_func(xR(1:7)); % obtain the pre-TD swing(left) toe position
        z_sw_pre_td = x_sw_pre_td(3);
        x_sw_pre_td = x_sw_pre_td(1);
        x_st_td = x_sw_pre_td; % <<< update the stance foot position for next step
    
        % update DRS position
        x_drs_pre_td = DRSmotion2_h(t_ode(end,1),1,1);
        x_drs_td = x_drs_pre_td; % <<< update the position of drs when touchdown
    
        % update the robot state xR by switching the left and right leg states
        xR_pre_td = xR;
        % position, angle can just switch
        x_COM_pre_td = xR_pre_td(1:3); 
        x_right_pre_td = xR_pre_td(4:5); 
        x_left_pre_td = xR_pre_td(6:7);
        % velocity, angular velocity need the impact reset map
        [qR_dot_pos_td_before_switch, impact_force] = impact_map(xR_pre_td,DRSmotion2_h,t_sim_f);
        x_COM_dot_pre_td_before_switch   = qR_dot_pos_td_before_switch(1:3);
        x_right_dot_pre_td_before_switch = qR_dot_pos_td_before_switch(4:5);
        x_left_dot_pre_td_before_switch  = qR_dot_pos_td_before_switch(6:7);
        xR_pos_td = [x_COM_pre_td; 
                     x_left_pre_td; 
                     x_right_pre_td;
                     x_COM_dot_pre_td_before_switch; 
                     x_left_dot_pre_td_before_switch; 
                     x_right_dot_pre_td_before_switch];
        % enforce the post-impact right foot position is exactly 0 height
        xR_pos_td(2) = xR_pos_td(2) - z_sw_pre_td;

        % % check the post-impact right foot height is 0
        % z_stance_pos_td = pos_toe_r_func(xR_pos_td(1:7));
        % 
        % % check the post-impact right foot velocity is 0
        % qR_dot_pos_td = xR_pos_td(8:14);
        % Jc_pos_td = sJcb_toe_r_func(xR_pos_td(1:7));
        % right_toe_vel_td = Jc_pos_td * qR_dot_pos_td;
        % right_toe_xz_vel_td = right_toe_vel_td([1 3],1);
        % 
        % % enforce the post-impact right foot velocity is 0
        % DRS_vx = DRSmotion2_h(t_td,1,2);
        % DRS_vz = DRSmotion2_h(t_td,3,2);
        % DRS_vRy = DRSmotion2_h(t_td,5,2);
        % DRS_plus = [DRS_vx; DRS_vz; DRS_vRy];
        % q_plus = Jc_pos_td([1 3 5],:) \ DRS_plus;
        
        xR = xR_pos_td;
    end

end % end while

%% 
%%%%%%%%%%%%%%%%%%
%%% check plot %%%
%%%%%%%%%%%%%%%%%%

% figure; hold on;
% plot(t_list, x_list(1,:),'-','LineWidth',2);
% plot(t_list, x_list(2,:),'-','LineWidth',2);
% plot(t_list, x_list(3,:),'-','LineWidth',2);
% plot(t_list, x_list(4,:),'-','LineWidth',2);
% plot(t_list, x_list(5,:),'-','LineWidth',2);
% plot(t_list, x_list(6,:),'-','LineWidth',2);
% plot(t_list, x_list(7,:),'-','LineWidth',2);
% legend('x','z','\theta','q1r','q2r','q1l','q2l')

% % plot snapshot
figure;%('Position',[-1427.8000 -2.2000 560 420]);
number_of_snapshot = 10;
ani_index = 1:round(length(t_list)/(number_of_snapshot-1)):length(t_list);
% ani_index = [ani_index, length(t_list)];
for i = 1:length(ani_index)
    index_to_plot = ani_index(i);
    bn = 1 - (i-1)/number_of_snapshot;
    plot_FiveLinkWalker(x_list(1:7,index_to_plot),bn);
end
plot([-1 1],[0 0],'k-o')
grid on;

% % plot ALIP state and step length decision for each loop
figure%('Position',[-827.8000 -2.2000 560 420]);
subplot(3,1,1);
plot(t2_list2, qA_list2(1,:),'LineWidth',2);
ylabel('x_{sc} (m)')
ylim([-0.3 0.3])

subplot(3,1,2);
plot(t2_list2, qA_list2(2,:),'LineWidth',2);
ylabel('L_s (kg*m^2/s)')
ylim([-50 50])

subplot(3,1,3);
plot(t2_list2, step_length_list2,'-','LineWidth',2);
xlabel('t (sec)')
ylabel('desired step length (m)')
ylim([-0.5 0.5])

sgtitle('ALIP states and step length controller')

% % plot task space desired and tracking result
name_list = ["c_z (m)", "\theta_b (rad)", "x_{sw} (m)", "z_{sw} (m)"];
y_low = [0.40 0.1 -0.5 -0.1];
y_upp = [0.80 0.3 0.2   0.1];
figure%('Position',[873 185 560 420]);
for i = 1:4
    subplot(2,2,i); hold on;
    plot(t2_list2,hc_list2(i,:),'r-','LineWidth',1);
    plot(t2_list2,hd_list2(i,:),'b-','LineWidth',1);
    plot(t2_list2,hd_cur_list2(i,:),'k--','LineWidth',1);
    plot(t2_list2(1:end-1),diff(step_list2),'-k')
    xlabel('t (sec)')
    ylabel(name_list(i));
    ylim([y_low(i) y_upp(i)])
end
legend({'robot state','desired (no DRS)','desired (DRS)'},'Location','best')

% % plot torque output
figure;
plot(t_list, u_list,'-','LineWidth',2)
legend('tau 1 right','tau 2 right','tau 1 left','tau 2 left')
ylim([-500 500])

%% plot x, z, theta
% figure('Position',[899.4000 183.4000 537.6000 524.8000]);
% subplot(3,1,1);
% plot(t_list, x_list(1,:),'-','LineWidth',2);
% ylabel('x (m)')
% subplot(3,1,2);
% plot(t_list, x_list(2,:),'-','LineWidth',2);
% ylabel('z (m)')
% subplot(3,1,3);
% plot(t_list, x_list(3,:),'-','LineWidth',2);
% ylabel('\theta (rad)')
% xlabel('time (s)')

%% plot q1r, q2r, q1l, q2l
% figure;
% subplot(2,2,1)
% plot(t_list, x_list(6,:),'-','LineWidth',2);
% ylabel('q1 left (rad)')
% subplot(2,2,3)
% plot(t_list, x_list(7,:),'-','LineWidth',2);
% ylabel('q2 left (rad)')
% subplot(2,2,2)
% plot(t_list, x_list(4,:),'-','LineWidth',2);
% ylabel('q1 right (rad)')
% subplot(2,2,4)
% plot(t_list, x_list(5,:),'-','LineWidth',2);
% ylabel('q2 right (rad)')
% 
% 
% figure;
% subplot(1,2,1); hold on;
% plot(t_list, x_list(6,:),'o-','LineWidth',2);
% plot(t_list, x_list(4,:),'s-','LineWidth',2);
% legend('q1 left (rad)','q1 right (rad)')
% 
% subplot(1,2,2); hold on;
% plot(t_list, x_list(7,:),'o-','LineWidth',2);
% plot(t_list, x_list(5,:),'s-','LineWidth',2);
% legend('q2 left (rad)','q2 right (rad)')

%% plot q1r_dot, q2r_dot, q1l_dot, q2l_dot
% figure;
% subplot(2,2,1)
% plot(t_list, x_list(6+7,:),'-','LineWidth',2);
% ylabel('q1 left dot (rad)')
% subplot(2,2,3)
% plot(t_list, x_list(7+7,:),'-','LineWidth',2);
% ylabel('q2 left dot (rad)')
% subplot(2,2,2)
% plot(t_list, x_list(4+7,:),'-','LineWidth',2);
% ylabel('q1 right dot (rad)')
% subplot(2,2,4)
% plot(t_list, x_list(5+7,:),'-','LineWidth',2);
% ylabel('q2 right dot (rad)')
% 
% figure;
% subplot(1,2,1); hold on;
% plot(t_list, x_list(6+7,:),'o-','LineWidth',2);
% plot(t_list, x_list(4+7,:),'s-','LineWidth',2);
% legend('q1 dot left (rad)','q1 dot right (rad)')
% subplot(1,2,2); hold on;
% plot(t_list, x_list(7+7,:),'o-','LineWidth',2);
% plot(t_list, x_list(5+7,:),'s-','LineWidth',2);
% legend('q2 dot left (rad)','q2 dot right (rad)')

%%
%%%%%%%%%%%%%%%%%
%%% animation %%%
%%%%%%%%%%%%%%%%%
close all
clear F
save = 0 ; mov=0; aviname='Five-link robot on DRS with ALIP controller 2'; % 設定要不要輸出avi檔案

figure('Position',[371.4000 109 1.1288e+03 583.2000]);
for i = 1:40:length(t2_list2) %size(x_traj,2)
    clf; 

    % animation
    subplot(3,4,[1 2 3 5 6 7 9 10 11])
    plot_FiveLinkWalker(x_list2(1:7,i));
    x_DRS = DRSmotion2_h(t2_list2(i),1,1);
    plot( -1+x_DRS:0.5:5+x_DRS ,(-1:0.5:5)*0,'k-o') % DRS
    grid on;
    axis equal
    xlim([-0.5 4.0])
    ylim([-1.2 2.2])
    % show info
    text(0,2,    ['t = ',num2str(t2_list2(i),'%.2f'),' sec'])
    text(0,1.75, ['DRS x = ',num2str(x_DRS,'%.5f'),' m'])
    text(2,2,    ['Desired velocity = ',num2str(controller.v_d,'%.2f'),' m/sec'])

    % ALIP state info: x_sc
    subplot(3,4,4); hold on;
    plot(t2_list2, qA_list2(1,:),'LineWidth',2);
    plot(t2_list2(i), qA_list2(1,i),'ro','LineWidth',1)
    ylabel('x_{sc} (m)')
    ylim([-0.3 0.3])
    title('ALIP states and step length controller')
    grid on;

    % ALIP state info: L_s
    subplot(3,4,8); hold on;
    plot([t2_list2(1) t2_list2(end)], [1 1]*controller.L_s_d,'r--','LineWidth',2);
    plot(t2_list2, qA_list2(2,:),'LineWidth',2);
    plot(t2_list2(i), qA_list2(2,i),'ro','LineWidth',1);
    ylabel('L_s (kg*m^2/s)')
    ylim([-20 20])
    grid on;

    % ALIP controller: step length
    subplot(3,4,12); hold on;
    plot(t2_list2, step_length_list2,'LineWidth',2);
    plot(t2_list2(i), step_length_list2(i),'ro','LineWidth',1);
    xlabel('t (sec)')
    ylabel('desired step length (m)')
    ylim([-0.5 0.5])
    grid on;    
    
    drawnow
    if save==1; mov=mov+1; F(mov)=getframe(1); end
end

if save==1
    writerObj = VideoWriter([aviname,'.avi']);
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for j=1:length(F)
        % convert the image to a frame
        frame = F(j) ;    
        if isempty(frame.cdata)
            continue;
        end
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end



%% subfunction

function qA = Robot2ALIP_state(xR)
% convert the xR (Full order robot states) to the qA (reduced-order ALIP
% model's states). 
%
% INPUT:
%   xR = [qR; qR_dot] [14x1]
%       qR = [x, z, theta, q1_r, q2_r, q1_l, q2_l]' [7x1]
%       (x, z, and theta are w.r.t. World frame)
%       qR_dot = d(qR)/dt [7x1]
%
% OUTPUT:
%   qA = [x_sc; L_s] [2x1]
%       x_sc: horizontal distance from stance toe position to CoM
%       position w.r.t. World frame {W}
%       L_s : Angular momentum about the stance toe position w.r.t.
%       Wrold frame {W}

qR = xR(1:7,1);
qR_dot = xR(8:14,1);
ps = pos_toe_r_func(qR); % position of the stance toe

CoM = pos_CoM_func(qR); % function generated by gen_five_link_model_usingFROST.m
x_sc = CoM(1) - ps(1);

AM = AMworld_about_ps_func(qR,qR_dot,ps); % function generated by gen_five_link_model_usingFROST.m
L_s = AM(2);

qA = [x_sc; L_s];

end

function hc = Robot2task_state(xR)
% convert the xR (Full order robot states) to the hc (task space states). 
% 
% INPUT:
%   xR = [qR; qR_dot] [14x1]
%       qR = [x, z, theta, q1_r, q2_r, q1_l, q2_l]' [7x1]
%       (x, z, and theta are w.r.t. World frame)
%       qR_dot = d(qR)/dt [7x1]
%
% OUTPUT:
%   hc = [z_COM, theta_b, x_sw, z_sw]'
%       swing foot = left foot
qR = xR(1:7,1);

% CoM = pos_CoM_func(qR);
% z_CoM = CoM(3);
% 
% theta_b = qR(3);
% 
% swing_toe = pos_toe_l_func(qR);
% x_sw = swing_toe(1);
% z_sw = swing_toe(3);
% 
% hc = [z_CoM; theta_b; x_sw; z_sw];
hc = hc_func(qR);
end

function [value,isterminal,direction] = SE_touchdown(t,x,DRSmotion_h)
% The Switch Event used in ode45 to terminate the integration once the
% swing leg (Left leg) touchdown
% 
% WE KEEP RIGHT FOOT AS SUPPORTING FOOT!
%         LEFT  FOOT AS SWING FOOT!
% 
% INPUT:
%   t: time
%   x: state [14x1]
%   DRSmotion_h: function handle of (t) and should return
%       DRS = [P,V,A] [6x3]
%           P: Position/Angle of the DRS, (m or rad), [6x1]
%               P = [Px, Py, Pz, Rx, Ry, Rz]'
%           V: Velocity of the DRS, (m/s or rad/s), [6x1]
%           A: Acceleration of the DRS, (m/s^2 or rad/s^2), [6x1]
%
% OUTPUT:
% value(i) is a mathematical expression describing the ith event. An event occurs when value(i) is equal to zero.
% isterminal(i) = 1 if the integration is to terminate when the ith event occurs. Otherwise, it is 0.
% direction(i) = 0 if all zeros are to be located (the default). A value of +1 locates only zeros where the event function is increasing, and -1 locates only zeros where the event function is decreasing. Specify direction = [] to use the default value of 0 for all events.
% check: https://www.mathworks.com/help/matlab/math/ode-event-location.html


q = x(1:7,1);
% q_dot = x(8:14,1);

left_toe_position = pos_toe_l_func(q);
% left_toe_Jacobian = sJcb_toe_l_func(q);

value = left_toe_position(3); % the height of the swing(left) toe
% should be modified if we consider vertical and pitching DRS motion

isterminal = 1;
direction = -1; % only trigger when height is decreasing

end

























