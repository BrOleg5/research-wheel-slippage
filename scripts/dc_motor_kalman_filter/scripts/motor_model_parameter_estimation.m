%% Initial parameter
R_a = 5.95; % Ohm
L_a = 8.9e-3; % H
k_t = 0.0514; % N*m/A
k_e = k_t;
J = 71e-7; % kg*m^2

V_nom = 24; % V
% I_nom = 0.86; % A
% omega_nom = 3600 * 2 * pi / 60; % rad/s

b = 1e-3; % N*m*s/rad
T_0 = 0.01; % N*m;

k_p = 2;

tau = 0.153; % s
%% Initial approximation
J = 1.4629e-06; % kg*m^2
b = 2.5841e-05; % N*m*s/rad
T_0 = 0.0082; % N*m;

k_p = 1.1958;
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

motor_reference = timeseries(T.(mstr + "setvel"), T.t);
motor_velocity = timeseries(T.(mstr + "vel"), T.t);
motor_current = timeseries(T.(mstr + "cur"), T.t);
%% Syncronization reference and output data
ref_step_idx = find(abs(diff(motor_reference.Data)) > 18);
vel_step_idx = find(abs(diff(motor_velocity.Data)) > 18);
assert(length(ref_step_idx) == length(vel_step_idx), "ref_step_idx and vel_step_idx have the same length.");
for i = 1:length(ref_step_idx)/2
    idx = ref_step_idx(2*i-1):vel_step_idx(2*i-1);
    motor_reference.Data(idx) = 0;

    idx = ref_step_idx(2*i):vel_step_idx(2*i);
    motor_reference.Data(idx) = motor_reference.Data(idx(1)-1);
end

tau = mean(T.t(vel_step_idx) - T.t(ref_step_idx));
%% Configure estimation
% link: https://uk.mathworks.com/help/sldo/ug/estimate-model-parameter-values-code.html
mdl = "dc_motor";
open_system(mdl);

% Create an experiment
Exp = sdo.Experiment(mdl);

SpeedSetPoint = Simulink.SimulationData.Signal;
SpeedSetPoint.Name = "Motor speed set-point";
SpeedSetPoint.BlockPath = mdl + "/speed set-point";
SpeedSetPoint.PortType = "outport";
SpeedSetPoint.PortIndex = 1;
SpeedSetPoint.Values = motor_reference;

Current = Simulink.SimulationData.Signal;
Current.Name = "Motor current";
Current.BlockPath = mdl + "/Current integrator";
Current.PortType = "outport";
Current.PortIndex = 1;
Current.Values = motor_current;

Velocity = Simulink.SimulationData.Signal;
Velocity.Name = "Motor velocity";
Velocity.BlockPath = mdl + "/Velocity integrator";
Velocity.PortType = "outport";
Velocity.PortIndex = 1;
Velocity.Values = motor_velocity;

Exp.InputData = SpeedSetPoint;
Exp.OutputData = [Current; Velocity];
%% Specify the Parameters to Estimate
% ParameterNames = ["tau", "k_p", "J", "b", "T_0", "L_a", "R_a", "k_t", "k_e"];
ParameterNames = ["k_p", "J", "b", "L_a", "R_a", "k_t", "k_e"];
ParameterIdx = 1:length(ParameterNames);
pd = containers.Map(ParameterNames, ParameterIdx);

p = sdo.getParameterFromModel(mdl, ParameterNames);

% tau
% d = 0.5;
% p(pd("tau")).Minimum = (1 - d) * tau;
% p(pd("tau")).Maximum = (1 + d) * tau;
% k_p
p(pd("k_p")).Minimum = 0;
p(pd("k_p")).Maximum = 10;
% J
p(pd("J")).Minimum = 0;
p(pd("J")).Maximum = 1e-2;
% b
p(pd("b")).Minimum = 0;
p(pd("b")).Maximum = 1e-1;
% T_0
% p(pd("T_0")).Minimum = 0;
% p(pd("T_0")).Maximum = 1e-2;
% L_a
d = 0.5;
p(pd("L_a")).Minimum = (1 - d) * L_a;
p(pd("L_a")).Maximum = (1 + d) * L_a;
% R_a
p(pd("R_a")).Minimum = (1 - d) * R_a;
p(pd("R_a")).Maximum = (1 + d) * R_a;
% k_t
p(pd("k_t")).Minimum = (1 - d) * k_t;
p(pd("k_t")).Maximum = (1 + d) * k_t;
% k_e
p(pd("k_e")).Minimum = (1 - d) * k_e;
p(pd("k_e")).Maximum = (1 + d) * k_e;

v = p;
%% Estimate the Parameters
Simulator = createSimulator(Exp);
Simulator = sim(Simulator);
% Define the Estimation Objective Function
estFcn = @(v) dc_motor_Objective(v, Simulator, Exp, pd);

set_param(mdl, 'MaxStep', '15e-3', 'FastRestart', 'on');
rng('default');
opt = sdo.OptimizeOptions('Method', 'lsqnonlin');
vOpt = sdo.optimize(estFcn, v, opt)
%% Compare the Measured Output and the Final Simulated Output
% Exp = setEstimatedValues(Exp, vOpt);
Exp = setEstimatedValues(Exp, v);

Simulator = createSimulator(Exp);
Simulator = sim(Simulator);
SimLog = find(Simulator.LoggedData, get_param(mdl, "SignalLoggingName"));
CurrentSignal = find(SimLog, "Motor current");
VelocitySignal = find(SimLog, "Motor velocity");
figure("Name", "Model simulation", "WindowState", "maximize");
subplot(2, 1, 1);
plot(motor_current, 'LineWidth', 1, 'Color', 'blue');
hold on;
plot(CurrentSignal.Values.Time, CurrentSignal.Values.Data, 'LineWidth', 1, 'Color', 'red');
grid on;
xlabel("Time, s");
ylabel("Current, A");
legend(["Measured", "Simulated"], 'FontSize', 14, 'Location', 'best');
title("Motor current");
ax = gca;
ax.FontSize = 14;
ax.XLimitMethod = 'tight';
subplot(2, 1, 2);
plot(motor_velocity, 'LineWidth', 1, 'Color', 'blue');
hold on;
plot(VelocitySignal.Values.Time, VelocitySignal.Values.Data, 'LineWidth', 1, 'Color', 'red');
grid on;
xlabel("Time, s");
ylabel("Velocity, rad/s");
legend(["Measured", "Simulated"], 'FontSize', 14, 'Location', 'best');
title("Motor velocity");
ax = gca;
ax.FontSize = 14;
ax.XLimitMethod = 'tight';
%% Validation
% Load experimental data
T = readtable("data\validation_data.csv", "NumHeaderLines", 2);

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

motor_reference = timeseries(T.(mstr + "setvel"), T.t);
motor_velocity = timeseries(T.(mstr + "vel"), T.t);
motor_current = timeseries(T.(mstr + "cur"), T.t);

tau = mean(T.t(vel_step_idx) - T.t(ref_step_idx));
% Configure estimation
% link: https://uk.mathworks.com/help/sldo/ug/estimate-model-parameter-values-code.html
mdl = "dc_motor";
open_system(mdl);

% Create an experiment
Exp = sdo.Experiment(mdl);

SpeedSetPoint = Simulink.SimulationData.Signal;
SpeedSetPoint.Name = "Motor speed set-point";
SpeedSetPoint.BlockPath = mdl + "/speed set-point";
SpeedSetPoint.PortType = "outport";
SpeedSetPoint.PortIndex = 1;
SpeedSetPoint.Values = motor_reference;

Current = Simulink.SimulationData.Signal;
Current.Name = "Motor current";
Current.BlockPath = mdl + "/Current integrator";
Current.PortType = "outport";
Current.PortIndex = 1;
Current.Values = motor_current;

Velocity = Simulink.SimulationData.Signal;
Velocity.Name = "Motor velocity";
Velocity.BlockPath = mdl + "/Velocity integrator";
Velocity.PortType = "outport";
Velocity.PortIndex = 1;
Velocity.Values = motor_velocity;

Exp.InputData = SpeedSetPoint;
Exp.OutputData = [Current; Velocity];

% Compare the Measured Output and the Final Simulated Output
Simulator = createSimulator(Exp);
Simulator = sim(Simulator);
SimLog = find(Simulator.LoggedData, get_param(mdl, "SignalLoggingName"));
CurrentSignal = find(SimLog, "Motor current");
VelocitySignal = find(SimLog, "Motor velocity");
figure("Name", "Model simulation");
subplot(2, 1, 1);
plot(motor_current, CurrentSignal.Values.Time, CurrentSignal.Values.Data);
grid on;
xlabel("Time, s");
ylabel("Current, A");
legend(["Measured", "Simulated"]);
title("Motor current");
subplot(2, 1, 2);
plot(motor_velocity, VelocitySignal.Values.Time, VelocitySignal.Values.Data);
grid on;
xlabel("Time, s");
ylabel("Velocity, rad/s");
legend(["Measured", "Simulated"]);
title("Motor velocity");
%% Local functions
function vals = dc_motor_Objective(v, Simulator, Exp, paramDict)
% DC_MOTOR_OBJECTIVE
%
%    The DC_MOTOR_OBJECTIVE function is used to compare model
%    outputs against experimental data.
%
%    vals = DC_MOTOR_OBJECTIVE(v, Simulator, Exp) 
%
%    The |v| input argument is a vector of estimated model parameter values
%    and initial states.
%
%    The |Simulator| input argument is a simulation object used 
%    simulate the model with the estimated parameter values.
%
%    The |Exp| input argument contains the estimation experiment data.
%
%    The |vals| return argument contains information about how well the
%    model simulation results match the experimental data and is used by
%    the |sdo.optimize| function to estimate the model parameters.
%
%    https://uk.mathworks.com/help/sldo/ug/estimate-model-parameters-per-experiment-code.html
%
%    See also sdo.optimize, sdoExampleCostFunction

arguments
    v param.Continuous,
    Simulator sdo.SimulationTest,
    Exp sdo.Experiment,
    paramDict containers.Map
end

r = sdo.requirements.SignalTracking('Method', 'Residuals', 'Normalize', 'on');

% Update the experiments with the estimated parameter values.
Exp  = setEstimatedValues(Exp, v);

% Simulate the model and compare model outputs with measured experiment data.
Simulator = createSimulator(Exp, Simulator);
Simulator = sim(Simulator);
SimLog  = find(Simulator.LoggedData, get_param(Exp.ModelName, 'SignalLoggingName'));
Error = [];
for ctSig=1:numel(Exp.OutputData)
    Signal = find(SimLog, Exp.OutputData(ctSig).Name);

    err = evalRequirement(r, Signal.Values, Exp.OutputData(ctSig).Values);
    
    Error = [Error; err(:)];
end

% Anti-oscillation constraint
% L_a = v(paramDict("L_a")).Value;
% R_a = v(paramDict("R_a")).Value;
% b = v(paramDict("b")).Value;
% J = v(paramDict("J")).Value;
% k_t = v(paramDict("k_t")).Value;
% k_p = v(paramDict("k_p")).Value;
% D = (L_a .* b + J .* R_a).^2 - 4 .* L_a .* J .* (R_a .* b + k_t.^2 + k_t .* k_p);
% Cleq = exp(D);

% Return the residual errors to the optimization solver.
vals.F = Error(:);
% vals.Cleq = Cleq(:);
end