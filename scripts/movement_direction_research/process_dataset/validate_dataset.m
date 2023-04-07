function res = validate_dataset(T)
% VALIDATE_DATASET Validate table with dataset.
%   res = VALIDATE_DATASET(T);

arguments
    T (:, 25) table
end
res = true;

assert(all(T.t >= 0), "Time must be positive");

% Test increasing time variable
positive_dt = diff(T.t) > 0;
if(any(~positive_dt))
    warning(strcat("Time must be increasing. Wrong rows:", sprintf(" %d", find(~positive_dt))));
    res = false;
end

% Test speed set point vector amplitude
treshold = 1e-6; % m/s
set_speed_vec = sqrt(T.xsetspeed.^2 + T.ysetspeed.^2);
is_right_speed = diff(set_speed_vec) < treshold;
set_speed_vec_mode = round(mode(set_speed_vec), 1);
if(any(~is_right_speed))
    warning(strcat(sprintf("Speed vector must be equal %f m/s. ", set_speed_vec_mode), ...
                           "Wrong rows:", sprintf(" %d", find(~is_right_speed))));
    res = false;
end

% Test zero rotational velocity of robot
is_zero_velocity = T.setvel == 0;
if(any(~is_zero_velocity))
    warning(strcat("Robot velosity must be equal 0 rad/s. Wrong rows: ", ...
            sprintf(" %d", find(~is_zero_velocity))));
    res = false;
end

% Test variables of motors
for m = 1:3
    mstr = "m" + num2str(m);

    % Test motor velocity set point
    mot_set_vel_limit = 130; % rad/s
    is_less_set_vel = abs(T.(mstr + "setvel")) < mot_set_vel_limit;
    if(any(~is_less_set_vel))
           warning(strcat(sprintf("Absolute value of motor %d set velosity must be less than %f rad/s. ", ...
                          m, mot_set_vel_limit), "Wrong rows:", sprintf(" %d", find(~is_less_set_vel))));
           res = false;
    end

    % Test motor actual velocity
    is_less_vel = abs(T.(mstr + "vel")) < mot_set_vel_limit;
    if(any(~is_less_vel))
           warning(strcat(sprintf("Absolute value of motor %d velosity must be less than %f rad/s. ", ...
                          m, mot_set_vel_limit), "Wrong rows:", sprintf(" %d", find(~is_less_vel))));
           res = false;
    end

    % Test motor actual current
    mot_cur_limit = 3.3; % A
    is_less_cur = abs(T.(mstr + "cur")) < mot_cur_limit;
    if(any(~is_less_cur))
           warning(strcat(sprintf("Absolute value of motor %d current must be less than %f A. ", ...
                          m, mot_cur_limit), "Wrong rows:", sprintf(" %d", find(~is_less_cur))));
           res = false;
    end
end

% Test robot global position
% X-axis
low_x_limit = 200; % pixels
high_x_limit = 1400; % pixels
in_limits = (T.xpos > low_x_limit) & (T.xpos < high_x_limit);
if(any(~in_limits))
    warning(strcat(sprintf("X axis position must be in interval [%f; %f] pixels. ", ...
                           low_x_limit, high_x_limit), ...
                   "Wrong rows:", sprintf(" %d", find(~in_limits))));
    res = false;
end
% Y-axis
low_y_limit = 100; % pixels
high_y_limit = 1000; % pixels
in_limits = (T.ypos > low_y_limit) & (T.ypos < high_y_limit);
if(any(~in_limits))
    warning(strcat(sprintf("Y axis position must be in interval [%f; %f] pixels. ", ...
                           low_y_limit, high_y_limit), ...
                   "Wrong rows:", sprintf(" %d", find(~in_limits))));
    res = false;
end
% Angle
in_limits = (T.ang >= 0) & (T.ang <= 360);
if(any(~in_limits))
     warning(strcat("Robot angle must be in interval [0; 360]. Wrong rows:", ...
                    sprintf(" %d", find(~in_limits))));
     res = false;
end
end