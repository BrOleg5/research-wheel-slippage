function K = calc_linearization_coefficient(varargin)
% CALC_LINEARIZATION_COEFFICIENT  Calculate linearization coefficients for
% relative wheel slippage.
%
%    K = CALC_LINEARIZATION_COEFFICIENT(x_axis_robot_speed, y_axis_robot_speed, robotKinematics)
%    Calculate linearization coefficients for relative wheel slippages
%    from X and Y axes speed in m/s (x_axis_robot_speed, y_axis_robot_speed).
%    robotKinematics is OMNITHREEWHEELEDROBOTKINEMATICS object.
%
%    K = CALC_LINEARIZATION_COEFFICIENT(movement_angle, robotKinematics)
%    Calculate linearization coefficients for relative wheel slippages from movement direction in rad.
%    robotKinematics is OMNITHREEWHEELEDROBOTKINEMATICS object.

switch nargin
    case 2
        mustBeNumeric(varargin{1});
        mustBeVector(varargin{1});
        assert(isa(varargin{2}, "OmniThreeWheeledRobotKinematics"), ...
            "Last argument must be OmniThreeWheeledRobotKinematics object.");
        K = impl_calc_linearization_coefficient(varargin{1}, varargin{2});
    case 3
        mustBeNumeric(varargin{1});
        mustBeVector(varargin{1});
        mustBeNumeric(varargin{2});
        mustBeVector(varargin{2});
        assert(isa(varargin{3}, "OmniThreeWheeledRobotKinematics"), ...
            "Last argument must be OmniThreeWheeledRobotKinematics object.");
        assert(length(varargin{1}) == length(varargin{2}), ...
            "First and second arguments must have the same length");
        angle = atan2(varargin{2}, varargin{1});
        K = impl_calc_linearization_coefficient(angle, varargin{3});
end
end
%% Local functions
function K = impl_calc_linearization_coefficient(movement_angle, robotKinematics)
arguments
    movement_angle {mustBeNumeric},
    robotKinematics (1, 1) OmniThreeWheeledRobotKinematics
end

movement_angle = anyvec2rowvec(movement_angle);
V = [cos(movement_angle); sin(movement_angle); zeros(size(movement_angle))];
M = robotKinematics.M * robotKinematics.wheel_radius;

n = length(movement_angle);
K = zeros(3, 1, n);
for i = 1:n
    K(:, :, i) = abs(M * V(:, i));
end
end