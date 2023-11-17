function x_dot = FiveLink_constrained_dynamics(t,x,u,DRSmotion)
% Compute the change of state variables of the five link robot with on foot
% constrained on the moving surface
%
% WE KEEP RIGHT FOOT AS SUPPORTING FOOT!
%         LEFT  FOOT AS SWING FOOT!
% 
% EOM:
%               D(q)*q_ddot + C_vec = B*u
% Constrain:
%               J_foot*q_ddot + J_foot_dot*q_dot ( = x_foot_ddot ) = As 
%
% INPUT:
%   t: time (sec) [1x1]
%   x: state, x = [q', q_dot']', [14x1]
%           q = [x, z, theta, q1_r, q2_r, q1_l, q2_l]' , [7x1]
%           q_dot = [...], [7x1]
%   u: control input, [4x1]
%       B, [7x4], constant input-selection matrix
%       u = [tau_q1_r, tau_q2_r, tau_q1_l, tau_q2_l] [4x1]
%   DRSmotion: function handle of (t) and should return
%       DRS = [P,V,A] [6x3]
%           P: Position/Angle of the DRS, (m or rad), [6x1]
%               P = [Px, Py, Pz, Rx, Ry, Rz]'
%           V: Velocity of the DRS, (m/s or rad/s), [6x1]
%           A: Acceleration of the DRS, (m/s^2 or rad/s^2), [6x1]
%
% OUTPUT:
%   x_dot: time derivative of the state variables

qR=x(1:7,1);
qR_dot=x(8:14,1);

% Call the original Equation of motion 
[D,C_vec,B] = FiveLink_dynamics(x);

% Call the Jacobian from the q to supporting point/foot
J_foot = sJcb_toe_r_func(qR); % [6x7] = [ (6DoF in 3D space) x (7DoF of Robot)]
J_foot_dot = sJcb_toe_r_dot_func(qR,qR_dot); 

% Call the DRS motion
DRS = DRSmotion(t);
a_DRS = DRS(1:6,3);

%%% reduce the a_DRS, Jacobian to a planar ones
J_foot = J_foot([1 3],:);
J_foot_dot = J_foot_dot([1 3],:);
a_DRS = a_DRS([1 3],:);

% Construct the C_vec_bar, B_bar
C_vec_bar = (C_vec)-J_foot'/(J_foot/D*J_foot')*(J_foot/D*(C_vec)-J_foot_dot*qR_dot+a_DRS);
B_bar = B-J_foot'/(J_foot/D*J_foot')*(J_foot/D*B);

q_ddot = D\(B_bar*u - C_vec_bar);

x_dot = [qR_dot; q_ddot];

end





