%% Load experimental data
T = readtable("data\data_for_identification_without_integral.csv", "NumHeaderLines", 2);

% Process data
T.t = T.t / 1000;
T.t = T.t + 0.03;
start_row = zeros(1, width(T));
T = [array2table(start_row, "VariableNames", T.Properties.VariableNames); T];

for m = 1:3
    mstr = "m" + num2str(m);
    T.(mstr + "cur") = sign(T.(mstr + "vel")) .* T.(mstr + "cur");
end

m = 1;
mstr = "m" + num2str(m);

time = T.t;
motor_reference = T.(mstr + "setvel");
motor_velocity = T.(mstr + "vel");
motor_current = T.(mstr + "cur");

%% Syncronization reference and output data
ref_step_idx = find(abs(diff(motor_reference)) > 18);
vel_step_idx = find(abs(diff(motor_velocity)) > 18);
assert(length(ref_step_idx) == length(vel_step_idx), "ref_step_idx and vel_step_idx have the same length.");
time_delay = mean(time(vel_step_idx) - time(ref_step_idx));
for i = 1:length(ref_step_idx)/2
    idx = ref_step_idx(2*i-1):vel_step_idx(2*i-1);
    motor_reference(idx) = 0;

    idx = ref_step_idx(2*i):vel_step_idx(2*i);
    motor_reference(idx) = motor_reference(idx(1)-1);
end
%% Inretpolate (resample) data with fixed sample rate
% link: https://uk.mathworks.com/help/matlab/ref/interp1.html
Ts = 1e-3;
desiredT = 0:Ts:time(end);
% Reference signal
interpolated_reference = interp1(time, motor_reference, desiredT);
% plot(time, motor_reference, '.-', desiredT, interpolated_reference, '.-');
% Current
interpolated_current = interp1(time, motor_current, desiredT);
% plot(time, motor_current, '.-', desiredT, interpolated_current, '.-');
% legend('Original', 'Resampled');
% Velocity
interpolated_velocity = interp1(time, motor_velocity, desiredT);
% plot(time, motor_velocity, '.-', desiredT, interpolated_velocity, '.-');
% legend('Original', 'Resampled');

interpolated_torque = -sign(interpolated_reference);

interpolated_reference = interpolated_reference.';
interpolated_torque = interpolated_torque.';
interpolated_current = interpolated_current.';
interpolated_velocity = interpolated_velocity.';
%% Create identification data
data = iddata([interpolated_current, interpolated_velocity], [interpolated_reference, interpolated_torque], ...
              Ts, 'ExperimentName', 'Motor identification', 'Name', 'Experimental data', ...
             'InputName', {'velocity set-point', 'torque'}, 'InputUnit', {'rad/s', 'N*m'}, ...
             'OutputName', {'current', 'velocity'}, 'OutputUnit', {'A', 'rad/s'});

% idplot(data);
%% Initial parameter
R_a = 5.95; % Ohm
L_a = 8.9e-3; % H
k_t = 0.0514; % N*m/A
k_e = k_t;
J = 71e-7; % kg*m^2

% V_nom = 24; % V
% I_nom = 0.86; % A
% omega_nom = 3600 * 2 * pi / 60; % rad/s

b = 1e-3; % N*m*s/rad
T_0 = 0.01; % N*m;

k_p = 2;
%% Configure Estimable Parameter of Grey-Box Model
odefun = @motor_state_space_model;

parameters = {'armature inductive', L_a; 'armature resistance', R_a; 'total moment of inertia', J;...
              'coefficient of viscous friction', b; 'torque constant', k_t; 'counter EMF constant', k_e; ...
              'speed P-controller coefficient', k_p; 'torque', T_0};

init_sys = idgrey(odefun, parameters, 'c');

ParameterNames = ["L_a", "R_a", "J", "b", "k_t", "k_e", "k_p", "T_0"];
ParameterIdx = 1:length(ParameterNames);
pd = containers.Map(ParameterNames, ParameterIdx);

% init_sys.Structure.Parameters(pd("L_a")).Free = false;
% init_sys.Structure.Parameters(pd("R_a")).Free = false;
% init_sys.Structure.Parameters(pd("k_t")).Free = false;
% init_sys.Structure.Parameters(pd("k_e")).Free = false;

init_sys.Structure.Parameters(pd("J")).Minimum = 0;
init_sys.Structure.Parameters(pd("b")).Minimum = 0;
init_sys.Structure.Parameters(pd("k_p")).Minimum = 0;
init_sys.Structure.Parameters(pd("T_0")).Minimum = 0;

d = 0.05;
init_sys.Structure.Parameters(pd("L_a")).Minimum = (1 - d) * init_sys.Structure.Parameters(pd("L_a")).Value;
init_sys.Structure.Parameters(pd("L_a")).Maximum = (1 + d) * init_sys.Structure.Parameters(pd("L_a")).Value;
init_sys.Structure.Parameters(pd("R_a")).Minimum = (1 - d) * init_sys.Structure.Parameters(pd("R_a")).Value;
init_sys.Structure.Parameters(pd("R_a")).Maximum = (1 + d) * init_sys.Structure.Parameters(pd("R_a")).Value;
init_sys.Structure.Parameters(pd("k_t")).Minimum = (1 - d) * init_sys.Structure.Parameters(pd("k_t")).Value;
init_sys.Structure.Parameters(pd("k_t")).Maximum = (1 + d) * init_sys.Structure.Parameters(pd("k_t")).Value;
init_sys.Structure.Parameters(pd("k_e")).Minimum = (1 - d) * init_sys.Structure.Parameters(pd("k_e")).Value;
init_sys.Structure.Parameters(pd("k_e")).Maximum = (1 + d) * init_sys.Structure.Parameters(pd("k_e")).Value;
%% Estimate
opt = greyestOptions('InitialState', 'zero', 'Display', 'on', 'DisturbanceModel', 'none', ...
                     'EnforceStability', true);
[est_sys, x_0] = greyest(data, init_sys, opt);
est_sys.Name = 'Estimated model';
%% Compare
opt = compareOptions('InitialCondition', 'zero');
compare(data, est_sys, Inf, opt);
grid on;
%% Simulate system
R_a = est_sys.Report.Parameters.ParVector(pd("R_a"));
L_a = est_sys.Report.Parameters.ParVector(pd("L_a"));
k_t = est_sys.Report.Parameters.ParVector(pd("k_t"));
k_e = est_sys.Report.Parameters.ParVector(pd("k_e"));
J = est_sys.Report.Parameters.ParVector(pd("J"));

b = est_sys.Report.Parameters.ParVector(pd("b"));
T_0 = est_sys.Report.Parameters.ParVector(pd("T_0"));

k_p = est_sys.Report.Parameters.ParVector(pd("k_p"));

% [A, B, C, D] = motor_state_space_model(L_a, R_a, J, b, k_t, k_e, k_p, T_0);
% sys = ss(A, B, C, D);
% lsim(sys, [motor_reference, motor_torque], time);
%% Local functions
function [A, B, C, D] = motor_state_space_model(L_a, R_a, J, b, k_t, k_e, k_p, T_0, Ts)
arguments
    L_a {mustBeNumeric},
    R_a {mustBeNumeric},
    J {mustBeNumeric},
    b {mustBeNumeric},
    k_t {mustBeNumeric},
    k_e {mustBeNumeric},
    k_p {mustBeNumeric},
    T_0 {mustBeNumeric},
    Ts {mustBeNumeric} = 0
end
if((L_a.*b + J.*R_a).^2 - 4.*L_a.*J.*(R_a.*b + k_t.^2 + k_t.*k_p) > 0)
    A = [-R_a/L_a, -(k_p + k_e)/L_a; 
            k_t/J,             -b/J];
    B = [k_p/L_a,     0;
               0, T_0/J];
else
    A = [-R_a/L_a, -(k_p + k_e)/L_a; 
            k_t/J,             -b/J];
    B = [k_p/L_a,     0;
               0, T_0/J];
end
C = [1, 0;
     0, 1];
D = zeros(2, 2);
end