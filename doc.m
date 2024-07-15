%% rotate_3D
%
% Function to perform rotation of a given vector or array of vectors around any given axis, in 2D or 3D spaces.
%
% Author : nicolas.douillet9 (at) gmail.com, 2017-2024.
%
%% Syntax
% R = rotate_3D(V, mode, theta);
%
% R = rotate_3D(V, mode, theta, u); 
%
% R = rotate_3D(V, mode, theta, u, angle_unit);
%
% [R,Rm] = rotate_3D(V, mode, theta, u, angle_unit);
%
%% Description
% R = rotate_3D(V, mode, theta) computes the vector R, which results
% from the rotation of V vector around one of the the basis vectors, which
% is choosen in the mode : 'x', 'y', 'z', 'i', 'j', 'k', or 'any'.
%
% R = rotate_3D(V, mode, theta, u) computes the vector R, which results
% from the rotation of V vector around u vector and of theta angle in radian.
%
% R = rotate_3D(V, mode, theta, u, angle_unit) uses angle_unit for theta
% unit (radian or degree).
%
% [R,Rm] = rotate_3D(V, mode, theta, u, angle_unit) also returns the
% rotation matrix.
%
% Important NB : in 2D -(xOy) plan- mandatory rotation axis is 'z' 'k'. It will
% be set as so by default if input is different. Also in 2D, in case u is missing it
% is automatically set to the origin [0,0]' by default.
%
%% Input arguments
%
%        [ -Vx- ]
% - V = [ -Vy- ], real (array of) vector(s) double, the vector(s) to rotate. Size(V) = [3,vector_nb]. 
%        [ -Vz- ]
%
% - mode : character string in the set {'x', 'X', 'y', 'Y', 'z', 'Z', 'i', 'I', 'j', 'J', 'k', 'K', 'any', 'ANY'}. 
%
% - theta : real scalar double, the rotation angle in radians or in signed degres.
%
%          [ux]
% - u : = [uy], real column vector double, the rotation axis.
%          [uz]
%
% - angle_unit : character string in the set {'radian','degree','RADIAN,'DEGREE'}. Case insensitive.
%
%% Output arguments
%
%        [ -Rx- ]
% - R = [ -Ry- ], real (array of) vector(s) double, the resulting rotated vector(s). Size(R) = [3,vector_nb].
%        [ -Rz- ]
%
% - Rm : real matrix double, the resulting rotation matrix, size(Rm) = [3,3].
% 
%% 3D Example #1
%
% 2*pi/3 rotation of basis vector i around [1 1 1]' vector is basis vector j.
V = [1 0 0]';
mode = 'any';
u = ones(3,1);
theta = 2*pi/3;
R = rotate_3D(V, mode, theta, u)
%
%% 3D Example #2
%
% 2*pi/3 rotation of basis vector i around [0 0 1]' / z axis. u argument
% is optional in this case.
mode = 'z';
R = rotate_3D(V, mode, theta) % equivalent to rotate_3D([1 0]', mode, theta) here
%
%% 3D Example #3
%
% Array example
% 2*pi/3 rotation of (i,j,k) basis vectors around [1 1 1]' is its permutation (j,k,i).
V = eye(3);
u = ones(3,1);
mode = 'any';
R = rotate_3D(V, mode, theta, u)
%
%% 3D Example #4
%
% Shifted rotation axis
% Regular hexagon in 3D space
V1 = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
V5 = [0 0 1]';
V2 = (2/3)*(V1+V3) - (1/3)*V5;
V4 = (2/3)*(V3+V5) - (1/3)*V1;
V6 = (2/3)*(V1+V5) - (1/3)*V3;
P = cat(2,V1,V2,V3,V4,V5,V6); % green hexagon
n = cross(V1-V5,V2-V5);
% 2*pi/3 rotation around line defined by the vector normal to P plane and going through V2
R = V2 + rotate_3D(P-V2, 'any', 2*pi/3, n); % red hexagon
figure;
line([P(1,:) P(1,1)],[P(2,:) P(2,1)],[P(3,:) P(3,1)],'Color',[0 1 0],'Linewidth',2), hold on;
line([R(1,:) R(1,1)],[R(2,:) R(2,1)],[R(3,:) R(3,1)],'Color',[1 0 0],'Linewidth',2), hold on;
line([V2(1) V2(1)+n(1)],[V2(2) V2(2)+n(2)],[V2(3) V2(3)+n(3)],'Color',[0 0 1],'Linewidth',2), hold on;
view(-45,27);
axis equal, axis tight;
%
%% 2D Example #1
%
% 3*pi/2 rotation of basis vector j around basis vector z is basis vector i.
V = [0 1]';
mode = 'z';
theta = 3*pi/2;
R = rotate_3D(V, mode, theta)
%
%% 2D Example #2
%
% pi/6 rotation of unitary vector [0.5 sqrt(3)/2]' around the origin is basis vector j.
% (default)
V = [0.5 sqrt(3)/2]';
mode = 'z';
theta = pi/6;
R = rotate_3D(V, mode, theta)
%
%% 2D Example #3
%
% 45° rotation of basis vector i around -0.5*[sqrt(2) sqrt(2)]' is basis vector j.
V = [1 0]';
u = -0.5*[sqrt(2) sqrt(2)]';
mode = 'x'; % 2D mode automatically set back to 'z' in case of mistake.
theta = 45; % unit i=in degree
R = rotate_3D(V, mode, theta, u, 'degree')