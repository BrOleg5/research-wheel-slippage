function K = calc_linearization_coefficient(x_axis_robot_speed, y_axis_robot_speed, robotKinematics)
% CALC_LINEARIZATION_COEFFICIENT

arguments
    x_axis_robot_speed {mustBeNumeric, mustBeVector},
    y_axis_robot_speed {mustBeNumeric, mustBeVector},
    robotKinematics OmniThreeWheeledRobotKinematics
end
n = length(x_axis_robot_speed);
assert(length(y_axis_robot_speed) == n, ...
       "x_axis_robot_speed and y_axis_robot_speed must have the same length");

angle = atan2(y_axis_robot_speed, x_axis_robot_speed);
omega = zeros(size(y_axis_robot_speed));
B = [cos(angle), sin(angle), sign(omega)];

M = robotKinematics.M * robotKinematics.wheel_radius;
K = zeros(3, 1, n);
for i = 1:n
    W = M * [x_axis_robot_speed(i); y_axis_robot_speed(i); omega(i)];
    K(:, :, i) = M * B(i, :).' .* sign(W);
end
end