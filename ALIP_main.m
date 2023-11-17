% ALIP simulation: using state-feedback stepping controller
% 2023/09/20

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Declear variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
m = 10; % kg
H = 0.5; % m
g = 9.81; % m/s^2

model_para.m = m;
model_para.H = H;
model_para.g = g;

% DRS motion
amp = 0.03; % m
freq = 3;  % rad/sec
Xs = @(t)amp*sin(freq*t); % DRS dynamics movement
Vs = @(t)amp*cos(freq*t)*freq; % DRS velocity

% controller setting, state-feedback version
T = 0.2; % sec, Step period
u_d = 0.1; % m, desired step length
x_d = [-0.08; 5]; % x_d = [x_sc_d; L_s_d]
K = [0.5 0.05]; % State feedback gain

controller.T = T;
controller.u_d = u_d;
controller.x_d = x_d;
controller.K = K;

%% plot Xs
t = 0:0.01:10;
Xs_t = Xs(t);
% figure; plot(t,Xs_t);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation setting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_sim = 10; % sec, total simulation length
x0 = [-0.1; % m, x_sc
      m*H*1];     % kg*m^2/s, L_s = m*Vc*H 

%%%%%%%%%%%%%%%%%%%%%%
%%% Run Simulation %%%
%%%%%%%%%%%%%%%%%%%%%%
x_traj = ALIP(x0, T_sim, Xs, Vs, model_para, controller);

%% plot

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

%% animation
close all

clear F
save=0; mov=0; aviname='ALIP'; % 設定要不要輸出avi檔案

figure('Position',[250 200 1500 800])
for i = 150:1:350 %size(x_traj,2)
    clf; 
    
    subplot(2,3,[1 2 3]);
    hold on;

    plot([-0.5 2.5]+Xs(x_traj(1,i)),[0 0],'k-','LineWidth',1.5) % ground

    plot(x_traj(4,i), H, 'ro','MarkerSize',20,'MarkerFaceColor','r')
    plot([x_traj(4:5,i)], [H 0],'k-','Linewidth',2);
    plot(x_traj(5,i), -0.01,'b^','LineWidth',2,'MarkerSize',5,'MarkerFaceColor','b')

    title(['ALIP t = ',num2str(x_traj(1,i),'%.3f')]);

    ylim([-0.05 1.0]);
    axis equal
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
    t_pre = t_ode(end,:) % check t_pre = tau_k = T_step
    x_pre = x_ode(end,:)';
    x_w_com = x_w_com0 + (x_ode(:,1)' - x_ode(1,1)) + Xs(t_ode') - Xs(t_ode(1));
    x_w_sup = ones(1,length(t_ode))*x_w_sup0        + Xs(t_ode') - Xs(t_ode(1));
    
    Dx = ALIP_foot_landing(x_pre,controller);
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

function Dx = ALIP_foot_landing(x_pre,controller)
% input:
%   x_pre [2x1]: pre-impact state variables
%       x1 = x_sc : COM x position relative to contact point on DRS
%       x2 = L_s  : angular momentum about contact point on DRS
%   t [1x1]: time of pre-impact
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
% output:
%   Dx [2x1]: state variables' jump during foot landing
%       Dx1 = Dx_sc
%       Dx2 = DL_s

Dx1 = - controller.K * x_pre ...
      + controller.K * controller.x_d ...
      - controller.u_d;
Dx2 = 0; % impact invariant

Dx = [Dx1; Dx2];
end






















