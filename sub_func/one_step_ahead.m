
function step_length = one_step_ahead(qA,t_cur,t_ctd,DRSmotion2_h,Robot,controller)
% given the ALIP state and the known DRSmotion, calculate a step_length for
% current state such that the angular momentum at next step is desired
% 
% INPUT:
%   qA: current ALIP state = [x_sc; L_s] [2x1]
%       x_sc: horizontal distance from stance toe position to CoM
%       position w.r.t. World frame {W}
%       L_s : Angular momentum about the stance toe position w.r.t.
%       Wrold frame {W}
%   t_cur [1x1]: current time
%   t_ctd [1x1]: coming touchdown time
%   DRSmotion2_h = @(t,r,c)DRSmotion_h = out [1x1]
%       out = DRS(r,c)
%           DRS = [P,V,A] [6x3]
%               P: Position/Angle of the DRS, (m or rad), [6x1]
%                   P = [Px, Py, Pz, Rx, Ry, Rz]'
%               V: Velocity of the DRS, (m/s or rad/s), [6x1]
%               A: Acceleration of the DRS, (m/s^2 or rad/s^2), [6x1]
%   Robot [structure]
%       .m : mass (kg)
%       .H : COM height (m)
%       .g : gravity (m/s^2)
%   controller [structure]
%       .L_s_d [1x1]: desired angular momentum, m*H*v_desired
%       .T [1x1]: duration for each step
% 
% OUTPUT:
%   step_length [1x1]

T_step = controller.T;
l = sqrt(Robot.g / Robot.H );
mHl = Robot.m * Robot.H * l;
delta_t = t_ctd - t_cur;

% compute L_s_km: AM about point s (supporting point) at T_k^-
Ax = [0, 1/(Robot.m*Robot.H);
      Robot.m*Robot.g, 0];

Vs_x = @(t)DRSmotion2_h(t,1,2); % DRS x velocity w.r.t. world frame
fx = @(t)[-Vs_x(t); 0];
exp_Ax_deltaT_fx = @(tau)expm(Ax*(t_ctd-tau))*fx(tau);
forcing_integral = integral(exp_Ax_deltaT_fx,t_cur,t_ctd,'ArrayValued',true);

T_t = [cosh(l*delta_t), sinh(l*delta_t)/mHl;
       mHl*sinh(l*delta_t), cosh(l*delta_t)];

x_t_km = T_t * qA + forcing_integral;
x_s_t_km = x_t_km(1);
L_s_km = x_t_km(2);

% compute Vx2_T_step:
exp_Ax_deltaT_fx_2 = @(tau)expm(Ax*(t_ctd+T_step-tau))*fx(tau);
forcing_integral_2 = integral(exp_Ax_deltaT_fx_2,t_ctd,t_ctd+T_step,'ArrayValued',true);
Vx2_T_step = forcing_integral_2(2);

% compute the position vector from "Swing foot position at T_km"
% to "COM position at T_km"
x_swc = (controller.L_s_d - Vx2_T_step - cosh(l*T_step)*L_s_km) / (mHl * sinh(l*T_step));

step_length = x_s_t_km - x_swc;

end