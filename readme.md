# Note
link to the hackMD
https://hackmd.io/@ichia1999/five_link_walker

# ALIP model

TBD

# Five-link robot

## Install FROST
Important link: 
https://github.com/ayonga/frost-dev
https://ayonga.github.io/frost-dev/index.html

Additional steps to install:
1. install Mathematica version<13.1, I install 13.0 and it works
2. copy additional four files to 
`frost-dev-master\third\mathlink`
3. follow the installation steps:
https://ayonga.github.io/frost-dev/pages/installation.html
(install: winGW, Mathematica Symbolic Toolbox for MATLAB--Version 2.0

Some degubing:
https://github.com/UMich-BipedLab/Cassie_Model/issues/5#issuecomment-1274028519
https://blog.csdn.net/qq_43309940/article/details/127040054 


## URDF file
Below are some useful reference for understanding the content of URDF file:
1. Tree structure, joint definition, only kinematics: 
http://wiki.ros.org/urdf/Tutorials/Create%20your%20own%20urdf%20file 
2. Dynamics (link property): 
http://wiki.ros.org/urdf/Tutorials/Adding%20Physical%20and%20Collision%20Properties%20to%20a%20URDF%20Model 
3. URDF/XML spec: 
http://wiki.ros.org/urdf/XML 
4. URDF Link properties:
Important to know the **COM position, inertia orientation**, ...
https://wiki.ros.org/urdf/XML/link
![](https://hackmd.io/_uploads/H1fKtKHfT.png)
5. URDF Joint properties:
Important to know the **length of each link, orientation of joint axis**, ...
https://wiki.ros.org/urdf/XML/joint
![](https://hackmd.io/_uploads/HJTSKYrM6.png)

## Dynamic and Kinematic model
Using FROST and URDF file to generate Kinematics and Equation of Motion (EOM) of five link robot.

### Goal of modeling
Given a URDF file, we need to
1. Understand the specfication of the robot, kinematic tree, the definition of the state $q$, the transformation between each link.
2. Obtain the forward kinematic function. That is the function to express the Cartesian position of each joint, end-effector, COM, and point of interest as a function of joint angle.
3. Obtain the Inertia matrix, Coriolis term, and Gravity vector to form the Equation of Motion (EOM) of the whole system. That is the $D(q), C(q,\dot{q})\dot{q}$, and $G(q)$ in the equation:
$D(q)\ddot{q}+C(q,\dot{q})\dot{q}+G(q)=\Gamma$

### Goal 1: Understand the URDF file content
First, include the FROST environment
```matlab=
clear
clsoe all
clc

% Step: 1
%  run >> frost_addpath
%         @TRACE Lab\research\project\digit_control\FROST_2\frost-dev-master
% addpath C:\Users\apple\OneDrive - purdue.edu\TRACE Lab\research\project\digit_control\FROST\frost-dev-master
% frost_addpath
addpath ..\FROST\frost-dev-master
frost_addpath
% After the window pop up
% --> Choose the C:\Program Files\Wolfram Research\Mathematica\13.0\math.exe file
```
Second, import the URDF file
```matlab=
%% Read .urdf file
FiveLink = RobotLinks('urdf\five_link_walker.urdf','planar');
% attribute:
%   planar: 2D for five link robot
%   floating: 6DoF for Digit
```
Show important information in the URDF file, which can help building the kinematic tree and the mapping between each link and joint's frame. For detail information, please refer to the previous section "URDF file"
```matlab= 
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
```

The result is shown in the below. 
```
Joint info.
--
Joint # 1 is BasePosX
From Parent link [Origin] --> [prismatic] joint --> Child link [BasePosX]
Joint Axis w.r.t. Parent's frame = 1  0  0
Offset (Child link's origin w.r.t. Parent's frame): 0  0  0
R (Child frame's orientation w.r.t. Parent's frame orientation
1  0  0
0  1  0
0  0  1
--
Joint # 2 is BasePosZ
From Parent link [BasePosX] --> [prismatic] joint --> Child link [BasePosZ]
Joint Axis w.r.t. Parent's frame = 0  0  1
Offset (Child link's origin w.r.t. Parent's frame): 0  0  0
R (Child frame's orientation w.r.t. Parent's frame orientation
1  0  0
0  1  0
0  0  1
--
Joint # 3 is BaseRotY
From Parent link [BasePosZ] --> [revolute] joint --> Child link [torso]
Joint Axis w.r.t. Parent's frame = 0  1  0
Offset (Child link's origin w.r.t. Parent's frame): 0  0  0
R (Child frame's orientation w.r.t. Parent's frame orientation
1  0  0
0  1  0
0  0  1
--
Joint # 4 is q1_right
From Parent link [torso] --> [continuous] joint --> Child link [right_thigh]
Joint Axis w.r.t. Parent's frame = 0  1  0
Offset (Child link's origin w.r.t. Parent's frame): 0  0  0
R (Child frame's orientation w.r.t. Parent's frame orientation
1  0  0
0  1  0
0  0  1
--
Joint # 5 is q2_right
From Parent link [right_thigh] --> [continuous] joint --> Child link [right_shin]
Joint Axis w.r.t. Parent's frame = 0  1  0
Offset (Child link's origin w.r.t. Parent's frame): 0           0         0.4
R (Child frame's orientation w.r.t. Parent's frame orientation
1  0  0
0  1  0
0  0  1
--
Joint # 6 is q1_left
From Parent link [torso] --> [continuous] joint --> Child link [left_thigh]
Joint Axis w.r.t. Parent's frame = 0  1  0
Offset (Child link's origin w.r.t. Parent's frame): 0  0  0
R (Child frame's orientation w.r.t. Parent's frame orientation
1  0  0
0  1  0
0  0  1
--
Joint # 7 is q2_left
From Parent link [left_thigh] --> [continuous] joint --> Child link [left_shin]
Joint Axis w.r.t. Parent's frame = 0  1  0
Offset (Child link's origin w.r.t. Parent's frame): 0           0         0.4
R (Child frame's orientation w.r.t. Parent's frame orientation
1  0  0
0  1  0
0  0  1
```
and
```
Link info.
--
Link # 1 is torso
Mass: 28
Inertia in the diagonal: 0        1.33           0
Inertia frame origin (COM position) w.r.t. Link frame: 0           0        0.24
Inertia frame orientation w.r.t. Link frame: 
1  0  0
0  1  0
0  0  1
--
Link # 2 is right_thigh
Mass: 1
Inertia in the diagonal: 0        0.47           0
Inertia frame origin (COM position) w.r.t. Link frame: 0           0        0.11
Inertia frame orientation w.r.t. Link frame: 
1  0  0
0  1  0
0  0  1
--
Link # 3 is left_thigh
Mass: 1
Inertia in the diagonal: 0        0.47           0
Inertia frame origin (COM position) w.r.t. Link frame: 0           0        0.11
Inertia frame orientation w.r.t. Link frame: 
1  0  0
0  1  0
0  0  1
--
Link # 4 is right_shin
Mass: 1
Inertia in the diagonal: 0         0.2           0
Inertia frame origin (COM position) w.r.t. Link frame: 0           0        0.24
Inertia frame orientation w.r.t. Link frame: 
1  0  0
0  1  0
0  0  1
--
Link # 5 is left_shin
Mass: 1
Inertia in the diagonal: 0         0.2           0
Inertia frame origin (COM position) w.r.t. Link frame: 0           0        0.24
Inertia frame orientation w.r.t. Link frame: 
1  0  0
0  1  0
0  0  1
```
These results can be used to check the definition of $q$ and the kinematic tree. That is 
$q = [x, z, \theta, q_{1r}, q_{2r}, q_{1l}, q_{2l}]$. The kinematic tree can be represented by the following figure:
![](https://hackmd.io/_uploads/rybZQcSG6.png)

To further understand the transformation between each link and the COM position of each link, the offset and R of link and joint are needed (also are shown above) and the result is expressed in the following figure:
![](https://lh3.googleusercontent.com/pw/ADCreHeMQ4ztBXHLEo4t-fatbEtraSJRdM4QO8TzoxHMH0sAtpYgb-Qgvwn2R2-_w28YjNe2dI8kLP6V1wBXgoUezzPwGKJKnvedeOxI-NmbgGCyiopU1RX8jdANeFvZVjKzB5VlpKq8liwmr84RfIBFFRLHkOQGnqFurmvOxUaY4qWWQhcVg4nRHMRozeUM71i8THAVESXb0wwvXqNllhEdX5UbpTP9f6ZCbd8orcWs-OklBK1kGeYiJOGNPWCuOV76YHVBs_QPoD6e92DuX8cVVQA1Y6GfCbJHGqxTd9XOPqq1SwVFQEbB3jESTweH3DzUTglP9wrcHtadVXnv1njQEqE5SCyNPnO6Wg0nEJMHA84LVGptppHfcIxWGnJ7jrzNdBvdVPXKsTqbFvgjSBV-Sq77W3ChgbHOCl8xzWO36crN7NZccFoP5m0ASdXqpRsN5_dyF40NhzZ2lbt7Scd9Aogks3lgTSWDNa_Y7PKu4hp1D0RdqcaKl2rq1fMoixapyghKRL4sXoBPVwhIZaIhfldbw1KHNsg8ZoqDAVDFms7rQcM2tw9pF5RSpO6nLBt3O6plTyvVrFlLhCvOWgAGwYRTNTU-uSxrc6g9Hnt3ui7a3XrQqf_iCn4TY4JBcn0CS2SCH_5sPd1xsobL1RODToJraMVNzXZqfsfrPSsgD9OobPseLa2R5TyHmum4f5O0DeXkLIqNp8fyraBSzUZoBi794NkN95gDjKqxmcVmikoWfgJR-YogWVeCxvncW0shiElcLTnAYi_xHRMBjJ699bF1-1zCmfMeNzY5JBTO2vKLjTMiD7OYE5EpwQyM1TnS6vvvnCsWwFU-IhyqdV0KJZmp0g7XXUtzx_VNuQ2rjDmHMSXuSV5UHReFIFy9nxk7mOkJPtt97Iwx32II2efhDApi8XM=w1323-h919-s-no-gm?authuser=0)

### Goal 2: Obtain the forward kinematics function

Obtain the COM position by using FROST function:
```matlab= !
%% Find COM position
pos_CoM = FiveLink.getComPosition();
% input: q [1x7]
% output: CoM position [1x3]
X = SymVariable('x',[7,1]); % five link: 3DoF on a plane, and 4 joint angle
pos_CoM_func = SymFunction('pos_CoM_func',pos_CoM,{X});
export(pos_CoM_func, 'gen_byFROST'); % save the function as a .mexw64 file at the folder 'gen_byFROST'
```

In addition, we can also get the knee joint position using:
```matlab=
X = SymVariable('x',[7,1]);

% Right knee joint
pos_knee_r = FiveLink.getCartesianPosition(FiveLink.Links(1,4).Reference);
pos_knee_r_func = SymFunction('pos_knee_r_func',pos_knee_r,{X});
export(pos_knee_r_func, 'gen_byFROST');
```
If the point of interest is not located at the origin of any given frame, like the position of head and toe, an offest $P$ w.r.t. the Link frame can be added to the getCartesianPosition function:
```matlab= !
% Right toe
P = [0,0,0.4]; % offset of the point of interest, here is the toe position w.r.t to the Frame of the right_shin
pos_toe_r = FiveLink.getCartesianPosition(FiveLink.Links(1,4).Reference,P);
pos_toe_r_func = SymFunction('pos_toe_r_func',pos_toe_r,{X});
export(pos_toe_r_func, 'gen_byFROST');
```
Since the Spacial Jacobian of the toe is important, we can also obtain the Jacobian using FROST:
```matlab=
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
```
The difference between the Body Jacobian and the Spacial Jacobian lies in the frame representing the velocity of the point of interest. For Body Jacobian, the velocity of point of interest is expressed in the body frame. For Spacial Jacobian, the velocity of the point of interest is expressed in the spacial frame, or world frame, or inertia frame, which is a more pouplar choice in the field of robotics. However, when study locomotion on DRS, expressing velocity using body frame might be a good choice as well.
ref: https://modernrobotics.northwestern.edu/nu-gm-book-resource/5-1-1-space-jacobian/#department


### Goal 3: Obtain the Equation of Motion (EOM)

To obtain the Inertia matrix, Coriolis term, and Gravity vector, i.e. $D(q), C(q,\dot{q})\dot{q}$, and $G(q)$, the following MATLAB code is used:

First, calculate the dynamics using FROST:

```matlab=
FiveLink.configureDynamics; % generate Dynamics
```
Then, export the function to folder:
```matlab=
%% Export Mass matrix
export(FiveLink.Mmat,'gen_byFROST')

%% Export C(q,dq)dq and G(q) vectors
for i = 1:22 % based on the size of FiveLink.Fvec
    export(FiveLink.Fvec{i},'gen_byFROST')
    % export C(q,dq)dq and G(q) vector to 'gen_byFROST' folder
    % C(q,dq)dq from Fvec{1:end-1, 1}
    % G(q)      from Fvec{end    , 1}
end

```

## Plot the robot
