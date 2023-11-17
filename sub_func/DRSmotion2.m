function out = DRSmotion2(t,DRSsetting,r,c)
% Compute the DRS motion given time t
% INPUT:
%   t: time (sec) [1x1]
%   r: row number
%   c: column number
% 
% OUTPUT:
%   out = DRS(r,c)
%       DRS = [P,V,A] [6x3]
%           P: Position/Angle of the DRS, (m or rad), [6x1]
%               P = [Px, Py, Pz, Rx, Ry, Rz]'
%           V: Velocity of the DRS, (m/s or rad/s), [6x1]
%           A: Acceleration of the DRS, (m/s^2 or rad/s^2), [6x1]

ampX = DRSsetting.ampX; freqX = DRSsetting.freqX;
ampY = DRSsetting.ampY; freqY = DRSsetting.freqY;
ampZ = DRSsetting.ampZ; freqZ = DRSsetting.freqZ;

ampRX = DRSsetting.ampRX; freqRX = DRSsetting.freqRX;
ampRY = DRSsetting.ampRY; freqRY = DRSsetting.freqRY;
ampRZ = DRSsetting.ampRZ; freqRZ = DRSsetting.freqRZ;

Px = ampX* sin(freqX*t);
Py = ampY* sin(freqY*t);
Pz = ampZ* sin(freqZ*t);
Rx = ampRX*sin(freqRX*t);
Ry = ampRY*sin(freqRY*t);
Rz = ampRZ*sin(freqRZ*t);

Vx =  ampX* freqX* cos(freqX*t);
Vy =  ampY* freqY* cos(freqY*t);
Vz =  ampZ* freqZ* cos(freqZ*t);
VRx = ampRX*freqRX*cos(freqRX*t);
VRy = ampRY*freqRY*cos(freqRY*t);
VRz = ampRZ*freqRZ*cos(freqRZ*t);

Ax  = -ampX * freqX*freqX* sin(freqX*t);
Ay  = -ampY * freqY*freqY* sin(freqY*t);
Az  = -ampZ * freqZ*freqZ* sin(freqZ*t);
ARx = -ampRX*freqRX*freqRX*sin(freqRX*t);
ARy = -ampRY*freqRY*freqRY*sin(freqRY*t);
ARz = -ampRZ*freqRZ*freqRZ*sin(freqRZ*t);

P = [Px Py Pz Rx Ry Rz]';
V = [Vx Vy Vz VRx VRy VRz]';
A = [Ax Ay Az ARx ARy ARz]';

DRS = [P, V, A];

out = DRS(r,c);
end
