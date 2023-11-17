function [qR_dot_pos_td_before_switch, impact_force] = impact_map(xR_pre_td,DRSmotion2_h,t_td)
% calculate the post impact velocity of robot states, given the pre-impact
% state and the DRS velocity.
% 
% INPUT:
%   xR_pre_td = [qR; qR_dot] [14x1]
%   DRSmotion2_h = @(t,r,c)DRSmotion_h = out [1x1]
%       out = DRS(r,c)
%           DRS = [P,V,A] [6x3]
%               P: Position/Angle of the DRS, (m or rad), [6x1]
%                   P = [Px, Py, Pz, Rx, Ry, Rz]'
%               V: Velocity of the DRS, (m/s or rad/s), [6x1]
%               A: Acceleration of the DRS, (m/s^2 or rad/s^2), [6x1]
%   t_td: time of touchdown [1x1]
%
% OUTPUT:
%   qR_dot_pos_td_before_switch: post-impact velocity,
%                            without change of leg index [7x1]
%   impact_force = [fx; fz] [2x1]: impact force during the impact event 

qR_pre_td = xR_pre_td(1:7);
qR_dot_pre_td = xR_pre_td(8:14);

[D,~,~] = FiveLink_dynamics(xR_pre_td);
Jc_all = sJcb_toe_l_func(qR_pre_td);
% only has constraints on x,z directions
Jc = Jc_all([1,3],:);

Mat = [D    , -Jc.';
       Jc   , zeros(2,2)];

DRS_vx = DRSmotion2_h(t_td,1,2);
DRS_vz = DRSmotion2_h(t_td,3,2);

Vec = [D*qR_dot_pre_td;
       DRS_vx;
       DRS_vz];

Vec_post = Mat \ Vec;

qR_dot_pos_td_before_switch = Vec_post(1:7);
impact_force = Vec_post(8:9);


end