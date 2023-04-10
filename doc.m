%% rotate_3D
%
% Function to compute the rotation of a vector in 2D and / or 3D spaces.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2023.
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
% is choosen in the mode : 'x', 'y', or 'z'.
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
% Important NB : in 2D -(xOy) plan- mandatory rotation axis is 'z'. It will
% be set as so by default if input is different. Also in 2D, in case u is missing it
% is automatically set to the origin [0,0]' by default.
%
%% Input arguments
%
%        [ -Vx- ]
% - V = [ -Vy- ], real (array of) vector(s) double, the vector(s) to rotate. Size(V) = [3,vector_nb]. 
%        [ -Vz- ]
%
% - mode : character string in the set {'x', 'X', 'y', 'Y', 'z', 'Z', 'any', 'ANY'}. 
%
% - theta : real scalar double, the rotation angle in radian.
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