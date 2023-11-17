function [hd,alpha_d] = desired_task_space_traj(step_length,x_sw_ini,x_st_ini,controller)
% calculate the desired task space reference trajectory (assuming no DRS
% motion, we consider the DRS motion latter)
%
% INPUT:
%   step_length: desired step length for this step
%   x_sw_ini: initial swing foot position w.r.t. world frame
%   x_st_ini: initial stance foot position w.r.t. world frame
%   controller [structure]
%       .BezierAlpha
%       .Bezier.x_sw_base
% 
% OUTPUT:
%   hd = @(s)[z_COM, theta_b, x_sw, z_sw]' [4x1]
%   alpha_d: desired alpha to generate the Bezier curve

alpha = controller.BezierAlpha;

x_sw_0 = x_sw_ini;
x_sw_6 = x_st_ini + step_length;
alpha_x_sw = x_sw_0 + controller.Bezier.x_sw_base * (x_sw_6 - x_sw_0);

alpha(:,3) = alpha_x_sw;

hd = @(s)Bezier_6th(alpha,s);
alpha_d = alpha;
end