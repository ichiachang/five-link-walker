function [h,dh_ds,ddh_dss]=Bezier_6th(alpha,s)
% Beizer polynomial
% 
% INPUT: 
%   alpha [7xn]: the coefficient
%       for each column = [a0; a1; a2; a3; a4; a5; a6]
%       a0 = h(s=0/6);
%       a1 = h(s=1/6);
%       a2 = h(s=2/6); ...
%       ...
%       a6 = h(s=6/6);
%   s [1xm]: phase variables
% 
% OUTPUT:
%   h [nx1]: position
%   dh_ds [nx1]: first derivative of the position (dh/ds)
%   ddh_dss [nx1]: second derivative of the position d(dh/ds)/ds
%   
    S=[
            (1-s).^6
            6*s.*(1-s).^5
            15*s.^2.*(1-s).^4
            20*s.^3.*(1-s).^3
            15*s.^4.*(1-s).^2
            6*s.^5.*(1-s).^1
            s.^6
        ]'; % [mx7]

    dS=[
            6*(s - 1).^5
            - 30*s.*(s - 1).^4 - 6.*(s - 1).^5
            30*s.*(s - 1).^4 + 60*s.^2.*(s - 1).^3
            - 60*s.^2.*(s - 1).^3 - 60*s.^3.*(s - 1).^2
            15*s.^4.*(2*s - 2) + 60*s.^3.*(s - 1).^2
            - 30*s.^4.*(s - 1) - 6*s.^5
            6*s.^5
        ]'; % [mx7]

    ddS=[
            30*(s - 1).^4
            - 120*s.*(s - 1).^3 - 60*(s - 1).^4
            240*s.*(s - 1).^3 + 30*(s - 1).^4 + 180*s.^2.*(s - 1).^2
            - 120*s.*(s - 1).^3 - 60*s.^3.*(2*s - 2) - 360*s.^2.*(s - 1).^2
            120*s.*(2*s - 2) + 180*s.^2.*(s - 1).^2 + 30*s.^4
            - 120*s.^3.*(s - 1) - 60*s.^4
            30*s.^4
        ]'; % [mx7]
    
      h      =   S*alpha; % [mx7] * [7xn] = [mxn]
     dh_ds   =  dS*alpha;
    ddh_dss  = ddS*alpha;

    h = h.';
    dh_ds = dh_ds.';
    ddh_dss = ddh_dss.';
end