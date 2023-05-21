%% Common
V_nom = 24; % V
I_nom = 0.86; % A
omega_nom = 3600 * 2 * pi / 60; % rad/s

R_a = 5.95; % Ohm
L_a = 8.9e-3; % H

k_t = 0.0514; % N*m/A
k_e = k_t;

J = 1.1407e-6; % kg*m^2

k_p = 1.6713;

tau = 0.153; % s
%% Motor 1
b = 2.3685e-5; % N*m*s/rad
T_0 = 0.0081312; % N*m
motor1 = struct('R_a', R_a, 'L_a', L_a, 'k_t', k_t, 'k_e', k_e, 'J', J, 'k_p', k_p, 'tau', tau, ...
                'b', b, 'T_0', T_0);
%% Motor 2
b = 2.0285e-5; % N*m*s/rad
T_0 = 0.0069712; % N*m
motor2 = struct('R_a', R_a, 'L_a', L_a, 'k_t', k_t, 'k_e', k_e, 'J', J, 'k_p', k_p, 'tau', tau, ...
                'b', b, 'T_0', T_0);
%% Motor 3
b = 2.3637e-5; % N*m*s/rad
T_0 = 0.0080917; % N*m
motor3 = struct('R_a', R_a, 'L_a', L_a, 'k_t', k_t, 'k_e', k_e, 'J', J, 'k_p', k_p, 'tau', tau, ...
                'b', b, 'T_0', T_0);