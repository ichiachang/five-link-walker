function [D,C_vec,B] = FiveLink_dynamics(x)
% Compute the Inertia Matrix, C vector (C(q,q_dot)q_dot), and B matrix in
% the Lagrange equation of a Five Link robot using the files generated by
% FROST. 
% The EOM is : D(q)*q_ddot + C_vec = B*u
%   where   q = [x, z, theta, q1_r, q2_r, q1_l, q2_l], [7x1]
%           D(q), [7x7], Inertia Matrix
%           C_vec = C(q,q_dot)*q_dot, [7x1]
%           B, [7x4], constant input-selection matrix
%           u = [tau_q1_r, tau_q2_r, tau_q1_l, tau_q2_l] [4x1]
% INPUT:
%       x = [q', q_dot']', state variables [14x1]
%       q = [x, z, theta, q1_r, q2_r, q1_l, q2_l]', [7x1]
%       q_dot = ..., [7x1]
% OUTPUT:
%       D(q), [7x7]
%       C_vec = C(q,q_dot)*q_dot, [7x1]
%       B, [7x1]

q=x(1:7,1);
dq=x(8:14,1);

% call the file generated by FROST
D = Mmat_five_link_walker(q);

c = zeros(7,1);    
c = Ce1_vec1_five_link_walker(q,dq)+c;
c = Ce1_vec2_five_link_walker(q,dq)+c;
c = Ce1_vec3_five_link_walker(q,dq)+c;
c = Ce1_vec4_five_link_walker(q,dq)+c;
c = Ce1_vec5_five_link_walker(q,dq)+c;
c = Ce1_vec6_five_link_walker(q,dq)+c;
c = Ce1_vec7_five_link_walker(q,dq)+c;

c = Ce2_vec1_five_link_walker(q,dq)+c;
c = Ce2_vec2_five_link_walker(q,dq)+c;
c = Ce2_vec3_five_link_walker(q,dq)+c;
c = Ce2_vec4_five_link_walker(q,dq)+c;
c = Ce2_vec5_five_link_walker(q,dq)+c;
c = Ce2_vec6_five_link_walker(q,dq)+c;
c = Ce2_vec7_five_link_walker(q,dq)+c;

c = Ce3_vec1_five_link_walker(q,dq)+c;
c = Ce3_vec2_five_link_walker(q,dq)+c;
c = Ce3_vec3_five_link_walker(q,dq)+c;
c = Ce3_vec4_five_link_walker(q,dq)+c;
c = Ce3_vec5_five_link_walker(q,dq)+c;
c = Ce3_vec6_five_link_walker(q,dq)+c;
c = Ce3_vec7_five_link_walker(q,dq)+c;

c = Ge_vec_five_link_walker(q,dq)+c;
C_vec = -c;  

B=zeros(7,4);
B(4:7,1:4)=eye(4);
end