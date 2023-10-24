% Generate the EOM of five-link-walking robot
% Date:2023/09/27

clear;
close all;
clc;

% Step: 1
%  run >> frost_addpath
%         @TRACE Lab\research\project\digit_control\FROST_2\frost-dev-master
% addpath C:\Users\apple\OneDrive - purdue.edu\TRACE Lab\research\project\digit_control\FROST\frost-dev-master
% frost_addpath
addpath ..\FROST\frost-dev-master
frost_addpath
% After the window pop up
% --> Choose the C:\Program Files\Wolfram Research\Mathematica\13.0\math.exe file

% Step: 2
%  run the following
%  ...

%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% Read .urdf file %%%
%%%%%%%%%%%%%%%%%%%%%%%
FiveLink = RobotLinks('urdf\five_link_walker.urdf','planar');
% attribute:
%   planar: 2D for five link robot
%   floating: 6DoF for Digit

%% show some important definition within the URDF file

% show the info of each joint
disp('Joint info.')
for i = 1:length(FiveLink.Joints)
    disp('--')
    disp(['Joint # ', num2str(i),' is ',FiveLink.Joints(1,i).Name])
    disp(['From Parent link [',FiveLink.Joints(1,i).Parent,'] --> [', ...
        FiveLink.Joints(1,i).Type, '] joint --> Child link [',FiveLink.Joints(1,i).Child,']']);
    disp(['Joint Axis w.r.t. Parent''s frame = ',num2str(FiveLink.Joints(1,i).Axis)]);
    disp(['Offset (Child link''s origin w.r.t. Parent''s frame): ',num2str(FiveLink.Joints(1,i).Offset)]);
    disp('R (Child frame''s orientation w.r.t. Parent''s frame orientation: ')
    disp(num2str(FiveLink.Joints(1,i).R))
end

% show the info of each link
disp('Link info.')
for i = 1:length(FiveLink.Links)
    disp('--')
    disp(['Link # ', num2str(i),' is ',FiveLink.Links(1,i).Name])
    disp(['Mass: ',num2str(FiveLink.Links(1,i).Mass)])
    diag_inertia = diag(FiveLink.Links(1,i).Inertia);
    disp(['Inertia in the diagonal: ',num2str(diag_inertia')])
%     disp(['Inertia: '])
%     disp(FiveLink.Links(1,i).Inertia)
    disp(['Inertia frame origin (COM position) w.r.t. Link frame: ',num2str(FiveLink.Links(1,i).Offset)])
    disp('Inertia frame orientation w.r.t. Link frame: ')
    disp(num2str(FiveLink.Links(1,i).R))
end

%%
%%%%%%%%%%%%%%%%%%
%%% Kinematics %%%
%%%%%%%%%%%%%%%%%%

%%% Find CoM position
pos_CoM = FiveLink.getComPosition();
% input: q [7x1]
% output: CoM position [1x3]
X = SymVariable('x',[7,1]); % five link: 3DoF on a plane, and 4 joint angle
pos_CoM_func = SymFunction('pos_CoM_func',pos_CoM,{X});
export(pos_CoM_func, 'gen_byFROST'); % save the function as a .mexw64 file at the folder 'gen_byFROST'

% test the generated function
%{
clear;
close all;
clc;

addpath gen_byFROST\
q = [1 1 0 0 0 0 0];
COM = pos_CoM_func(q) % you can call
%}

%% Find right/left knee and toe position, and head position
X = SymVariable('x',[7,1]);

% Right knee joint
pos_knee_r = FiveLink.getCartesianPosition(FiveLink.Links(1,4).Reference);
pos_knee_r_func = SymFunction('pos_knee_r_func',pos_knee_r,{X});
export(pos_knee_r_func, 'gen_byFROST');

% Left knee joint
pos_knee_l = FiveLink.getCartesianPosition(FiveLink.Links(1,5).Reference);
pos_knee_l_func = SymFunction('pos_knee_l_func',pos_knee_l,{X});
export(pos_knee_l_func, 'gen_byFROST');

% Right toe
P = [0,0,0.4]; % offset of the point of interest, here is the toe position w.r.t to the Frame of the right_shin
pos_toe_r = FiveLink.getCartesianPosition(FiveLink.Links(1,4).Reference,P);
pos_toe_r_func = SymFunction('pos_toe_r_func',pos_toe_r,{X});
export(pos_toe_r_func, 'gen_byFROST');

% Left toe
P = [0,0,0.4]; % offset of the point of interest, here is the toe position w.r.t to the Frame of the left_shin
pos_toe_l = FiveLink.getCartesianPosition(FiveLink.Links(1,5).Reference,P);
pos_toe_l_func = SymFunction('pos_toe_l_func',pos_toe_l,{X});
export(pos_toe_l_func, 'gen_byFROST');

% Head
P = [0,0,0.63];
pos_head = FiveLink.getCartesianPosition(FiveLink.Links(1,1).Reference,P);
pos_head_func = SymFunction('pos_head_func',pos_head,{X});
export(pos_head_func, 'gen_byFROST');

%% find Spacial Jacobian (w.r.t. world frame)
X = SymVariable('x',[7,1]);

% Right toe
P = [0,0,0.4];
sJcb_toe_r = FiveLink.getSpatialJacobian(FiveLink.Links(1,4).Reference,P);
sJcb_toe_r_func = SymFunction('sJcb_toe_r_func',sJcb_toe_r,{X});
export(sJcb_toe_r_func, 'gen_byFROST');

% Left toe
P = [0,0,0.4];
sJcb_toe_l = FiveLink.getSpatialJacobian(FiveLink.Links(1,5).Reference,P);
sJcb_toe_l_func = SymFunction('sJcb_toe_l_func',sJcb_toe_l,{X});
export(sJcb_toe_l_func, 'gen_byFROST');


%%
%%%%%%%%%%%%%%%%
%%% Dynamics %%%
%%%%%%%%%%%%%%%%
% Calculate Dynamic equations
FiveLink.configureDynamics; % generate Dynamics


%% Export Mass matrix
export(FiveLink.Mmat,'gen_byFROST')
% for i = 1:6 % based on the size of FiveLink.Mmat
%     export(FiveLink.Mmat{i},'gen_byFROST')
%     % export Mass matrix to 'gen_byFROST' folder
%     % M(q)      from Mmat
% end

%% Export C(q,dq)dq and G(q) vectors
for i = 1:22 % based on the size of FiveLink.Fvec
    export(FiveLink.Fvec{i},'gen_byFROST')
    % export C(q,dq)dq and G(q) vector to 'gen_byFROST' folder
    % C(q,dq)dq from Fvec{1:end-1, 1}
    % G(q)      from Fvec{end    , 1}
end



