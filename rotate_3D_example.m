% rotate_3D_example
%
% Author : nicolas.douillet9 (at) gmail.com, 2017-2024.


step = pi/36;
theta = (pi/2)*ones(181,1); % (pi/2)*ones(181,1); % (0:step:pi)';
phi = 0:step:2*pi;

Cx = cos(phi);
Cy = sin(phi);
Cz = zeros(size(phi));

C1 = [Cx; Cy; Cz];

C2 = rotate_3D(C1, 'x', pi/2);
C3 = rotate_3D(C1, 'y', -pi/2);

b = sqrt(3)/3*[1 1 1]'; % 1st bissectrice vector
S1 = zeros(3, 121);

% 1st third of the green circle
for k = 0:121
    S1(:,k+1) = rotate_3D([1 0 0]', 'any', 2*pi*k/363, b);
end

S2 = rotate_3D(S1, 'any', 2*pi/3, b); % 2nd third of the green circle
S3 = rotate_3D(S1, 'any', -2*pi/3, b); % 3rd third of the green circle

figure;
set(gcf, 'Color', [0 0 0]);
plot3(C1(1,:), C1(2,:), C1(3,:), 'Color', 'w', 'LineWidth', 2), hold on;
plot3(C2(1,:), C2(2,:), C2(3,:), 'Color', 'w', 'LineWidth', 2), hold on;
plot3(C3(1,:), C3(2,:), C3(3,:), 'Color', 'w', 'LineWidth', 2), hold on;

line([0 1], [0 0], [0 0], 'Color', 'r', 'LineWidth', 3), hold on;
line([0 0], [0 1], [0 0], 'Color', 'r', 'LineWidth', 3), hold on;
line([0 0], [0 0], [0 1], 'Color', 'r', 'LineWidth', 3), hold on;
line([0 b(1)], [0 b(2)], [0 b(3)], 'Color', 'g', 'LineWidth', 3), hold on;

plot3(S1(1,:), S1(2,:), S1(3,:), 'Color', 'g', 'LineWidth', 3), hold on;
plot3(S2(1,:), S2(2,:), S2(3,:), 'Color', 'g', 'LineWidth', 3), hold on;
plot3(S3(1,:), S3(2,:), S3(3,:), 'Color', 'g', 'LineWidth', 3), hold on;

set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
axis equal, axis square, axis tight;
view([102,17]);