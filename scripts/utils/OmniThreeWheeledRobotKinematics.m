classdef OmniThreeWheeledRobotKinematics
    % OMNITHREEWHEELEDROBOTKINEMATICS Represent kinematics of omnidirectional three-wheeled robot.
    % E.g. Festo Robotino
    
    properties
        wheel_radius (1, 1) {mustBeNumeric, mustBePositive} = 1, % Radius of omni-wheel in meters
        robot_radius (1, 1) {mustBeNumeric, mustBePositive} = 1, % Radius of robot in meters
        M (3, 3) {mustBeNumeric}, % Inverse kinematics matrix
        W (3, 3) {mustBeNumeric} % Direct kinematics matrix
    end
    
    methods
        function obj = OmniThreeWheeledRobotKinematics(wheel_radius, robot_radius)
            % OMNITHREEWHEELEDROBOTKINEMATICS Construct an instance of this class

            arguments
                wheel_radius (1, 1) {mustBeNumeric, mustBePositive},
                robot_radius (1, 1) {mustBeNumeric, mustBePositive}
            end

            obj.robot_radius = robot_radius;
            obj.wheel_radius = wheel_radius;
            obj.M = [  -sin(pi/3),   cos(pi/3), obj.robot_radius;
                         -sin(pi),     cos(pi), obj.robot_radius;
                     -sin(5*pi/3), cos(5*pi/3), obj.robot_radius] / obj.wheel_radius;
            obj.W = [      -2*sin(pi/3),         -2*sin(pi),     -2*sin(5*pi/3);
                            2*cos(pi/3),          2*cos(pi),      2*cos(5*pi/3);
                     1/obj.robot_radius, 1/obj.robot_radius, 1/obj.robot_radius] / 3 * obj.wheel_radius;
            assert(all((inv(obj.M) - obj.W) < 1e-5, 'all'), ...
                   "Inverse and direct kinematics matrix not consistent.");
        end
        
        function [wheel1_velocity, wheel2_velocity, wheel3_velocity] = ...
                inverse(obj, x_axis_robot_speed, y_axis_robot_speed, rotational_robot_speed)
            % INVERSE Calculate inverse kinematics of robot
            %
            %   [wheel1_velocity, wheel2_velocity, wheel3_velocity] = INVERSE(obj, x_axis_robot_speed, y_axis_robot_speed, rotational_robot_speed);
            %       x_axis_robot_speed and y_axis_robot_speed in m/s.
            %       rotational_robot_speed in rad/s.
            %       Method return array with velocity of three wheels in rad/s.
            %
            % See also DIRECT

            arguments
                obj OmniThreeWheeledRobotKinematics,
                x_axis_robot_speed {mustBeNumeric, mustBeVector},
                y_axis_robot_speed {mustBeNumeric, mustBeVector},
                rotational_robot_speed {mustBeNumeric, mustBeVector}
            end
            
            n = length(x_axis_robot_speed);
            assert(length(y_axis_robot_speed) == n, ...
                   "x_axis_robot_speed and y_axis_robot_speed must have same length.");
            assert(length(rotational_robot_speed) == n, ...
                   "x_axis_robot_speed and rotational_robot_speed must have same length.");
           
            sz = size(x_axis_robot_speed);
            wheel1_velocity = zeros(sz);
            wheel2_velocity = zeros(sz);
            wheel3_velocity = zeros(sz);
            for i = 1:n
                wheel_velocity = obj.M * [x_axis_robot_speed(i); y_axis_robot_speed(i); ...
                                          rotational_robot_speed(i)];
                wheel1_velocity(i) = wheel_velocity(1);
                wheel2_velocity(i) = wheel_velocity(2);
                wheel3_velocity(i) = wheel_velocity(3);
            end
        end

        function [x_axis_robot_speed, y_axis_robot_speed, rotational_robot_speed] = ...
            direct(obj, wheel1_velocity, wheel2_velocity, wheel3_velocity)
            % DIRECT Calculate direct kinematics of robot
            %
            %   [x_axis_robot_speed, y_axis_robot_speed, rotational_robot_speed] = DIRECT(obj, wheel1_velocity, wheel2_velocity, wheel3_velocity);
            %       Wheel velocity in rad/s.
            %       Method return array with robot speed.
            %       x_axis_robot_speed and y_axis_robot_speed in m/s,
            %       rotational_robot_speed in rad/s.
            %
            % See also INVERSE

            arguments
                obj OmniThreeWheeledRobotKinematics,
                wheel1_velocity {mustBeNumeric, mustBeVector},
                wheel2_velocity {mustBeNumeric, mustBeVector},
                wheel3_velocity {mustBeNumeric, mustBeVector}
            end
        
            n = length(wheel1_velocity);
            assert(length(wheel2_velocity) == n, ...
                   "wheel1_velocity and wheel2_velocity must have same length.");
            assert(length(wheel3_velocity) == n, ...
                   "wheel1_velocity and wheel3_velocity must have same length.");
           
            sz = size(wheel1_velocity);
            x_axis_robot_speed = zeros(sz);
            y_axis_robot_speed = zeros(sz);
            rotational_robot_speed = zeros(sz);
            for i = 1:n
                robot_speed = obj.W * [wheel1_velocity(i); wheel2_velocity(i); wheel3_velocity(i)];
                x_axis_robot_speed(i) = robot_speed(1);
                y_axis_robot_speed(i) = robot_speed(2);
                rotational_robot_speed(i) = robot_speed(3);
            end
        end
    end
end

