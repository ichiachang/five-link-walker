% ALIP simulation: using one-step ahead stepping controller
% 2023/11/9

clear
close all
% clc

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Declear variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
m = 32; % kg
H = 0.6372; % m
g = 9.81; % m/s^2

model_para.m = m;
model_para.H = H;
model_para.g = g;

% DRS motion
amp = 0.05; % m
period = 1; % sec
freq = 2*pi/period;  % rad/sec
Xs = @(t)amp*sin(freq*t); % DRS dynamics movement
Vs = @(t)amp*cos(freq*t)*freq; % DRS velocity

% controller setting, one-step-forward version
T = 0.35; % sec, Step period
v_des = 0.5;
L_s_d = m*H*v_des; % desired angular momentum, m*H*v_desired, let v_desired = 1 m/s

controller.T = T;
controller.L_s_d = L_s_d;

%% plot Xs
t = 0:0.01:10;
Xs_t = Xs(t);
% figure; plot(t,Xs_t);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation setting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_sim = 10; % sec, total simulation length
x0 = [0.08; % m, x_sc
      m*H*-0.90];     % kg*m^2/s, L_s = m*Vc*H 

%%%%%%%%%%%%%%%%%%%%%%
%%% Run Simulation %%%
%%%%%%%%%%%%%%%%%%%%%%
x_traj = ALIP(x0, T_sim, Xs, Vs, model_para, controller);

%% 
%%%%%%%%%%%%
%%% plot %%%
%%%%%%%%%%%%

figure('Position',[-1842 522 1342 420]); 
subplot(1,3,1); hold on;
plot(x_traj(1,:), x_traj(2,:),'o-');
xlabel('time (sec)');
ylabel('x_{sc} (m)');

subplot(1,3,2); hold on;
plot(x_traj(1,:), x_traj(3,:));
xlabel('time (sec)');
ylabel('L_{s} (kg*m^2/s)');

subplot(1,3,3); hold on;
plot(x_traj(2,:), x_traj(3,:));
xlabel('x_{sc} (m)');
ylabel('L_{s} (kg*m^2/s)');

figure;
plot(x_traj(1,:), x_traj(4,:),'o-');
xlabel('time')
ylabel('COM position w.r.t world frame')

figure;
plot(x_traj(1,:), x_traj(5,:),'o-');
xlabel('time')
ylabel('Foot position w.r.t world frame')

%%
%%%%%%%%%%%%%%%%%
%%% animation %%%
%%%%%%%%%%%%%%%%%

close all

clear F
save=0; mov=0; aviname='ALIP on DRS_One-step ahead_1109'; % 設定要不要輸出avi檔案

figure('Position',[250 200 1500 800])
for i = 1:5:250 %size(x_traj,2)
    clf; 
    
    subplot(2,3,[1 2 3]);
    hold on;

    plot([-0.5 2.5]+Xs(x_traj(1,i)),[0 0],'k-','LineWidth',1.5) % ground
    plot([-0.5:0.5:2.5]+Xs(x_traj(1,i)),zeros(1,7),'or')

    plot(x_traj(4,i), H, 'ro','MarkerSize',20,'MarkerFaceColor','r')
    plot([x_traj(4:5,i)], [H 0],'k-','Linewidth',2);
    plot(x_traj(5,i), -0.01,'b^','LineWidth',2,'MarkerSize',5,'MarkerFaceColor','b')

    title(['ALIP t = ',num2str(x_traj(1,i),'%.3f'), ...
        ', DRS amp = ',num2str(amp),' m, period = ',num2str(period),' sec']);    

    axis equal
    ylim([-0.05 1.0]);
    xlim([-1.5 3.5]);
    grid on;

    subplot(2,3,4); hold on;
    plot(x_traj(1,:), x_traj(2,:),'-');
    plot(x_traj(1,i), x_traj(2,i),'ro');
    xlabel('time (sec)');
    ylabel('x_{sc} (m)');
    
    subplot(2,3,5); hold on;
    plot(x_traj(1,:), x_traj(3,:));
    plot(x_traj(1,i), x_traj(3,i),'ro');
    xlabel('time (sec)');
    ylabel('L_{s} (kg*m^2/s)');
    
    subplot(2,3,6); hold on;
    plot(x_traj(2,:), x_traj(3,:));
    plot(x_traj(2,i), x_traj(3,i),'ro');
    xlabel('x_{sc} (m)');
    ylabel('L_{s} (kg*m^2/s)');
    
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

function x_traj = ALIP(x0,T_sim,Xs,Vs,model_para,controller)
% input:
%   x0 [2x1]: initial state variables
%       x0(1) = x_sc_0 
%       x0(2) = L_s_0
%   T_sim [1x1]: Total duration of simulation
%   Xs [1x1]: DRS motion (position) , Xs=@(t)...
%   Vs [1x1]: DRS motion (velocity) , Vs=@(t)...
%   model_para:
%       .m : mass (kg)
%       .H : COM height (m)
%       .g : gravity (m/s^2)
%   controller:
%       .T : step period
%       .u_d : desired step length
%       .x_d : desired pre-impact state
%       .K : state feedback gain
% output:
%   x_traj [2xn]: state trajectory of simulation


t_res = 0.005; % sec, resolution
x_traj = []; 

x_w_com0 = Xs(0) + x0(1);   % COM x position w.r.t. world frame
x_w_sup0 = 0;               % foot contact point w.r.t world frame
t_sim0 = 0;
T_step = controller.T;

while(t_sim0 < T_sim)
    % stance phase integration
    t_span = t_sim0:t_res:t_sim0+T_step;
    odefunc = @(t,x)ALIP_stance(x,t,Xs,Vs,model_para);
    [t_ode, x_ode] = ode45(odefunc, t_span, x0);

    % foot landing reset map
    t_pre = t_ode(end,:); % check t_pre = tau_k = T_step
    x_pre = x_ode(end,:)';
    x_w_com = x_w_com0 + (x_ode(:,1)' - x_ode(1,1)) + Xs(t_ode') - Xs(t_ode(1));
    x_w_sup = ones(1,length(t_ode))*x_w_sup0        + Xs(t_ode') - Xs(t_ode(1));
    
    %%% Compute the Step Length and determined the state jump of this step
    % Here, we use 'one-step-forward' controller
    % 
    % since this is perfect ALIP model, the result of one-step-forward
    % controller should be the same for the same stance step duration
    index_current_time = 50;
    x_cur = x_ode(index_current_time,:)';
    t_cur = t_ode(index_current_time,1);
    t_ctd = t_sim0 + T_step;
    Dx = ALIP_foot_landing(x_cur,t_cur,t_ctd,model_para,controller,Vs);

    x_post = x_pre + Dx; % post-impact state

    % generate the trajectory of this segment
    traj_k = [t_ode'; x_ode'; x_w_com; x_w_sup]; 
    x_traj = [x_traj, traj_k];

    % update variables for next while iteration
    t_sim0 = t_sim0 + T_step;
    x0 = x_post;
    x_w_com0 = x_w_com(end);
    x_w_sup0 = x_w_sup(end) - Dx(1);
end

end



function dx = ALIP_stance(x,t,Xs,Vs,model_para)
% input:
%   x [2x1]: state variables
%       x1 = x_sc : COM x position relative to contact point on DRS
%       x2 = L_s  : angular momentum about contact point on DRS
%   t [1x1]: time
%   Xs [1x1]: DRS motion (position) , Xs=@(t)...
%   Vs [1x1]: DRS motion (velocity) , Vs=@(t)...
%   model_para :
%       .m : mass (kg)
%       .H : COM height (m)
%       .g : gravity (m/s^2)
% output:
%   dx [2x1]: state variables' time derivate

dx1 = x(2)/(model_para.m * model_para.H) - Vs(t);
dx2 = model_para.m * model_para.g * x(1);

dx = [dx1; dx2];
end

function Dx = ALIP_foot_landing(x_cur,t_cur,t_ctd,model_para,controller,Vs)
% Use controller to decide the foot landing position
%   Here, we use one-step-forward controller
%   ref to Yuan's TRO_2023 paper
% 
% input:
%   x_cur [2x1]: current state variables
%       x1 = x_sc : COM x position relative to contact point on DRS
%       x2 = L_s  : angular momentum about contact point on DRS
%   t_cur [1x1]: current time
%   t_ctd [1x1]: coming touchdown time
%   Xs [1x1]: DRS motion, Xs=@(t)...
%   model_para :
%       .m : mass (kg)
%       .H : COM height (m)
%       .g : gravity (m/s^2)
%   controller:
%       .T : step period
%       .u_d : desired step length
%       .x_d : desired pre-impact state
%       .K : state feedback gain
%       .L_s_d: desired angular momentum, m*H*v_desired
%   Vs [1x1]: DRS motion (velocity) , Vs=@(t)...
% 
% output:
%   Dx [2x1]: state variables' jump during foot landing
%       Dx1 = Dx_sc
%       Dx2 = DL_s

T_step = controller.T;
l = sqrt(model_para.g / model_para.H );
mHl = model_para.m * model_para.H * l;
delta_t = t_ctd - t_cur;

% compute L_s_km: AM about point s (supporting point) at T_k^-
Ax = [0, 1/(model_para.m*model_para.H);
      model_para.m*model_para.g, 0];
fx = @(t)[-Vs(t); 0];
exp_Ax_deltaT_fx = @(tau)expm(Ax*(t_ctd-tau))*fx(tau);
forcing_integral = integral(exp_Ax_deltaT_fx,t_cur,t_ctd,'ArrayValued',true);

T_t = [cosh(l*delta_t), sinh(l*delta_t)/mHl;
       mHl*sinh(l*delta_t), cosh(l*delta_t)];

x_t_km = T_t * x_cur + forcing_integral;
x_s_t_km = x_t_km(1);
L_s_km = x_t_km(2);

% compute Vx2_T_step:
exp_Ax_deltaT_fx_2 = @(tau)expm(Ax*(t_ctd+T_step-tau))*fx(tau);
forcing_integral_2 = integral(exp_Ax_deltaT_fx_2,t_ctd,t_ctd+T_step,'ArrayValued',true);
Vx2_T_step = forcing_integral_2(2);

% compute the position vector from "Swing foot position at T_km"
% to "COM position at T_km"
x_swc = (controller.L_s_d - Vx2_T_step - cosh(l*T_step)*L_s_km) / (mHl * sinh(l*T_step));

u_step_length = x_s_t_km - x_swc;

Dx1 = - u_step_length; % step_length
Dx2 = 0; % impact invariant

Dx = [Dx1; Dx2];
end


% function out = forcing_integral(t_cur,t_ctd,model_para,Vs)
% % compute the integral of forcing term in the solution of ALIP dynamics
% % with DRS forcing input
% %
% % INPUT:
% %   t_cur [1x1]: current time
% %   t_ctd [1x1]: coming touchdown time
% %   model_para :
% %       .m : mass (kg)
% %       .H : COM height (m)
% %       .g : gravity (m/s^2)
% %   Vs [1x1]: DRS motion (velocity) , Vs=@(t)...
% % OUTPUT:
% %   out
% 
% Ax = [0, 1/(model_para.m*model_para.H);
%       model_para.m*model_para.g, 0];
% fx = [-Vs; 0];
% exp_Ax_deltaT_fx = @(t)exp(Ax*(t_ctd-t))*fx;
% 
% t_list = linspace(t_cur,t_ctd,100);
% to_be_trapzed = exp_Ax_deltaT_fx(t_list);
% out = trapz(t_list, to_be_trapzed);
% 
% end



















