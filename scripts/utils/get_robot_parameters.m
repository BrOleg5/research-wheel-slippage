function [robotParameters, robotKinematics] = get_robot_parameters()
% GET_ROBOT_PARAMETERS  Get struct with Robotino parameters and object of Robotino kinematics
% Resistance of motor winding
robotParameters.R_a = 5.95; % Ohm
% Inductance of motor winding
robotParameters.L_a = 8.9e-3; % H
% Torque constant
robotParameters.k_t = 0.0514;
% Сounter EMF equation constant
robotParameters.k_e = robotParameters.k_t;

robotParameters.wheel_radius = 40e-3; % m
% robotParameters.robot_radius = 130e-3; % m
robotParameters.robot_radius = 125e-3; % m
robotParameters.gear_ratio = 16;

robotKinematics = OmniThreeWheeledRobotKinematics(robotParameters.wheel_radius, ...
                                                  robotParameters.robot_radius);
end