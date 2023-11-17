function u = IO_linearization(hc,hd_cur,dhd_ds_cur,ddhd_dss_cur,xR,t_cur,DRSsetting,controller)
% input/output linearization low-level controller
% 
% INPUT:
%   hc: current task space state
%   hd_cur: current desired task space stete
%   dhd_ds_cur: current first derivative of desired task space stete to s
%   ddhd_dss_cur: current second derivative of desired task space stete
%                   w.r.t. s
%   xR: full order robot state [14x1]
%   t_cur: current time
%   DRSsetting [structure] 
%       .ampX/Y/Z/Rx/Ry/Rz
%       .freqX/Y/Z/Rx/Ry/Rz
%   controller [structure]
%       .T
%       .IO_kp
%       .IO_kd
% 
% OUTPUT:
%   u: control input for the state space model
%       = [tau_1r, tau_2r, tau_1l, tau_2l]

qR=xR(1:7,1);
qR_dot=xR(8:14,1);
s_dot = 1 / controller.T;

% Call the original Equation of motion 
[D,C_vec,B] = FiveLink_dynamics(xR);

% Call the Jacobian from the q to supporting point/foot
J_foot = sJcb_toe_r_func(qR); % [6x7] = [ (6DoF in 3D space) x (7DoF of Robot)]
J_foot_dot = sJcb_toe_r_dot_func(qR,qR_dot); % <<< need to be updated!!

% Call the DRS motion
DRS = DRSmotion(t_cur,DRSsetting);
a_DRS = DRS(1:6,3);

%%% reduce the a_DRS, Jacobian to a planar ones
J_foot = J_foot([1 3],:);
J_foot_dot = J_foot_dot([1 3],:);
a_DRS = a_DRS([1 3],:);

% Construct the C_vec_bar, B_bar
C_vec_bar = (C_vec)-J_foot'/(J_foot/D*J_foot')*(J_foot/D*(C_vec)-J_foot_dot*qR_dot+a_DRS);
B_bar = B-J_foot'/(J_foot/D*J_foot')*(J_foot/D*B);

% the Jacobian of output function sJcb_hc_func, sJcb_hc_dot_func
J_hc = sJcb_hc_func(qR);
J_hc_dot = sJcb_hc_dot_func(qR,qR_dot);

% output function y, y_dot
y = hc - hd_cur;
ds_dt = 1 / controller.T;
y_dot = J_hc*qR_dot - dhd_ds_cur * ds_dt;
Kp = [10 10 1 1].' * controller.IO_kp;
Kd = [10 10 1 1].' * controller.IO_kd;
v =  - Kp .* y - Kd .* y_dot;

% final output
u = (J_hc/(D)*B_bar)\(J_hc/(D)*C_vec_bar + ddhd_dss_cur*s_dot - J_hc_dot*qR_dot + v);

end