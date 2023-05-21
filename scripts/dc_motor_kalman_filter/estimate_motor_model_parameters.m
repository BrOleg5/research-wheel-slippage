function estimate_motor_model_parameters()
lang = "ru";
export_graph = true;
export_graph_folder = "output_files/graphs/dc_motor_kalman_filter";
export_graph_extensions = ["emf", "pdf"];

log_folder = "output_files/logs/dc_motor_kalman_filter";
if(~isfolder(log_folder))
    mkdir(log_folder);
end
log_file_name = fullfile(log_folder, "dc_motor_model_identification.txt");
logID = fopen(log_file_name, 'w');

mdl = "dc_motor";
open_system(mdl);

p = specify_estimated_parameters(mdl);

% Load dataset
train_data = readtable("datasets\dc_motor_model_identification\data_for_identification_without_integral.csv", ...
                       "NumHeaderLines", 2);
train_data = process_dataset(train_data);
[train_data, train_tau] = sync_reference_and_measurements(train_data);
fprintf(logID, "Train tau = %.8g s\n", train_tau);

test_data = readtable("datasets\dc_motor_model_identification\validation_data.csv", ...
                      "NumHeaderLines", 2);
test_data = process_dataset(test_data);
[test_data, test_tau] = sync_reference_and_measurements(test_data);
fprintf(logID, "Test tau = %.8g s\n", test_tau);

for motor_num = 1:3
    fprintf(logID, "Estimate model parameters of motor %d\n", motor_num);
    mstr = "m" + num2str(motor_num);
    
    % Parameter estimation
    train_motor_reference = timeseries(train_data.(mstr + "setvel"), train_data.t);
    train_motor_velocity = timeseries(train_data.(mstr + "vel"), train_data.t);
    train_motor_current = timeseries(train_data.(mstr + "cur"), train_data.t);
    train_exp = create_experiment(mdl, train_motor_reference, train_motor_current, ...
                                  train_motor_velocity);

    vOpt = estimate_parameters(mdl, train_exp, p, 'LogFileID', logID);

    model_parameters_path = "output_files/model_parameters";
    if(~isfolder(model_parameters_path))
        mkdir(model_parameters_path);
    end
    save(fullfile(model_parameters_path, "motor" + num2str(motor_num) + ".mat"), "vOpt");
    
    train_exp = setEstimatedValues(train_exp, vOpt);

    fprintf(logID, "\tTrain data:\n");
    train_fig = plot_compare_measured_and_simulated_output(mdl, train_exp, 'Language', lang, ...
                                                           'LogFileID', logID);
    if(export_graph)
        output_dir = fullfile(export_graph_folder, "motor_model_parameter_estimation", ...
                              "motor" + num2str(motor_num));
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, "train_data");
        export_graphs(train_fig, file_name, export_graph_extensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(train_fig);
    end

    % Validation
    test_motor_reference = timeseries(test_data.(mstr + "setvel"), test_data.t);
    test_motor_velocity = timeseries(test_data.(mstr + "vel"), test_data.t);
    test_motor_current = timeseries(test_data.(mstr + "cur"), test_data.t);
    test_exp = create_experiment(mdl, test_motor_reference, test_motor_current, ...
                                 test_motor_velocity);

    test_exp = setEstimatedValues(test_exp, vOpt);

    fprintf(logID, "\tTest data:\n");
    test_fig = plot_compare_measured_and_simulated_output(mdl, test_exp, 'Language', lang, ...
                                                          'LogFileID', logID);
    if(export_graph)
        file_name = fullfile(output_dir, "test_data");
        export_graphs(test_fig, file_name, export_graph_extensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(test_fig);
    end

    % Set initial parameters value of 2nd and 3rd motors as 1st motor
    % Optimize estimation time for 2nd and 3rd motors
    if(motor_num == 1)
        p = vOpt;
    end
end
fclose(logID);
end
%% Local functions
function out_T = process_dataset(in_T)
arguments
    in_T table
end

check_table_vars(in_T.Properties.VariableNames, ["t", "m1cur", "m2cur", "m3cur", ...
                                                 "m1vel", "m2vel", "m3vel"]);

out_T = in_T;
% Conver time from ms to s
out_T.t = out_T.t / 1000;

% Shift time by 0.03 s so that time start from 0
out_T.t = out_T.t + 0.03;

% Add start row to table. Fill this row with zeros
start_row = zeros(1, width(out_T));
out_T = [array2table(start_row, "VariableNames", out_T.Properties.VariableNames); out_T];

% Add sign to current values. Sign taken from motor/wheel velocity
for m = 1:3
    mstr = "m" + num2str(m);
    out_T.(mstr + "cur") = sign(out_T.(mstr + "vel")) .* out_T.(mstr + "cur");
end
end

function [out_T, tau] = sync_reference_and_measurements(in_T, options)
arguments
    in_T table,
    options.DeadZone = 3
end

check_table_vars(in_T.Properties.VariableNames, ["t", "m1setvel", "m2setvel", "m3setvel", ...
                                                 "m1vel", "m2vel", "m3vel"]);

out_T = in_T;
% Find indexes of reference step changing
ref_step_idx = find(diff(abs(out_T.m1setvel) < options.DeadZone));
% Find indexes of velocity step changing
vel_step_idx = find(diff(abs(out_T.m1vel) < options.DeadZone));

assert(length(ref_step_idx) == length(vel_step_idx), ...
       "ref_step_idx and vel_step_idx have the same length.");
assert(rem(length(ref_step_idx), 2) == 0, "ref_step_idx and vel_step_idx must have even length.");

for i = 1:length(ref_step_idx)/2
    rise_gap_idx = ref_step_idx(2*i-1):vel_step_idx(2*i-1);
    fall_gap_idx = ref_step_idx(2*i):vel_step_idx(2*i);
    for m = 1:3
        setvelstr = "m" + num2str(m) + "setvel";
        out_T.(setvelstr)(rise_gap_idx) = 0;
        out_T.(setvelstr)(fall_gap_idx) = out_T.(setvelstr)(fall_gap_idx(1)-1);
    end
end

tau = mean(out_T.t(vel_step_idx) - out_T.t(ref_step_idx));
end

function p = specify_estimated_parameters(mdl)
arguments
    mdl {mustBeTextScalar}
end
ParameterNames = ["k_p", "J", "b", "T_0", "L_a", "R_a", "k_t", "k_e"];
p = sdo.getParameterFromModel(mdl, ParameterNames);

% Create parameters map to access parameters from p by names
ParameterIdx = 1:length(ParameterNames);
pd = containers.Map(ParameterNames, ParameterIdx);

% Initial model parameters. Parameters taken from documtentation of DC motor Dunkermotoren GR 42x25
p(pd("R_a")).Value = 5.95;
p(pd("R_a")).Info.Unit = "Ohm";

p(pd("L_a")).Value = 8.9e-3;
p(pd("L_a")).Info.Unit = "H";

p(pd("k_t")).Value = 0.0514;
p(pd("k_t")).Info.Unit = "N*m/A";
% Assume k_e = k_t
p(pd("k_e")).Value = p(pd("k_t")).Value;
p(pd("k_e")).Info.Unit = "V*s/rad";

p(pd("J")).Value = 71e-7;
p(pd("J")).Info.Unit = "kg*m^2";

p(pd("b")).Value = 1e-3;
p(pd("b")).Info.Unit = "N*m*s/rad";

p(pd("T_0")).Value = 1e-3;
p(pd("T_0")).Info.Unit = "N*m";

p(pd("k_p")).Value = 2;
p(pd("k_p")).Info.Unit = "V*s/rad";

% Set tunable and fixed parameters
p(pd("R_a")).Free = false;
p(pd("L_a")).Free = false;
p(pd("k_t")).Free = false;
p(pd("k_e")).Free = false;
p(pd("J")).Free = true;
p(pd("b")).Free = true;
p(pd("T_0")).Free = true;
p(pd("k_p")).Free = true;

% Set parameter bounds
p(pd("k_p")).Minimum = 0;
p(pd("k_p")).Maximum = 10;
p(pd("J")).Minimum = 0;
p(pd("J")).Maximum = 1e-2;
p(pd("b")).Minimum = 0;
p(pd("b")).Maximum = 1e-1;
p(pd("T_0")).Minimum = 0;
p(pd("T_0")).Maximum = 1e-2;
end

function Exp = create_experiment(mdl, motor_reference, motor_current, motor_velocity)
arguments
    mdl {mustBeTextScalar},
    motor_reference timeseries,
    motor_current timeseries,
    motor_velocity timeseries
end

% link: https://uk.mathworks.com/help/sldo/ug/estimate-model-parameter-values-code.html
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
end

function vOpt = estimate_parameters(mdl, Exp, v, options)
arguments
    mdl {mustBeTextScalar},
    Exp sdo.Experiment,
    v param.Continuous,
    options.Method {mustBeTextScalar, mustBeMember(options.Method, ["lsqnonlin", "fmincon", ...
                                                                    "fminsearch", "patternsearch", ...
                                                                    "surrogateopt"])} = "lsqnonlin",
    options.LogFileID {mustBeNumeric, mustBeInteger} = -1
end
Simulator = createSimulator(Exp);
Simulator = sim(Simulator);
% Define the Estimation Objective Function
estFcn = @(v) dc_motor_Objective(v, Simulator, Exp);

set_param(mdl, 'MaxStep', '15e-3', 'FastRestart', 'on');
rng('default');
opt = sdo.OptimizeOptions('Method', options.Method);
[vOpt, optimInfo] = sdo.optimize(estFcn, v, opt);

if(options.LogFileID >= 0)
    fprintf(options.LogFileID, "\tOpimization/estimation method: %s\n\tAlgorithm: %s\n", options.Method, ...
            optimInfo.SolverOutput.output.algorithm);
    fprintf(options.LogFileID, "\tIterations: %d\n\tFunc-count: %d\n\tSquared norm of the residual: %f\n", ...
            optimInfo.iterations, optimInfo.Stats.FuncCount, optimInfo.SolverOutput.resnorm);
    fprintf(options.LogFileID, "\tEstimated parameters:\n");
    for i = 1:length(vOpt)
        if(vOpt(i).Free)
            status = "tunable";
        else
            status = "fixed";
        end
        fprintf(options.LogFileID, "\t\t%s = %.8g %s (%s)\n", vOpt(i).Name, vOpt(i).Value, ...
                vOpt(i).Info.Unit, status);
    end
end
end

function vals = dc_motor_Objective(v, Simulator, Exp)
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
    Exp sdo.Experiment
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
    Signal = getElement(SimLog, Exp.OutputData(ctSig).Name);

    err = evalRequirement(r, Signal.Values, Exp.OutputData(ctSig).Values);
    
    Error = [Error; err(:)];
end

% Return the residual errors to the optimization solver.
vals.F = Error(:);
end

function fig = plot_compare_measured_and_simulated_output(mdl, Exp, options)
arguments
    mdl {mustBeTextScalar},
    Exp sdo.Experiment,
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en",
    options.LogFileID {mustBeNumeric, mustBeInteger} = -1
end

Simulator = createSimulator(Exp);
Simulator = sim(Simulator);

SimLog = find(Simulator.LoggedData, get_param(mdl, "SignalLoggingName"));
simulated_current = getElement(SimLog, "Motor current");
simulated_velocity = getElement(SimLog, "Motor velocity");

exp_dataset = Simulink.SimulationData.Dataset;
for i = 1:numel(Exp.OutputData)
    exp_dataset = addElement(exp_dataset, Exp.OutputData(i));
end
measured_current = getElement(exp_dataset, "Motor current");
measured_velocity = getElement(exp_dataset, "Motor velocity");

resnorm = calculate_squared_norm_of_residual(SimLog, exp_dataset);
fprintf(options.LogFileID, "\t\tSquared norm of the residual: %g\n", resnorm);

fig = figure("Name", "Compare measured and simulated outputs", "WindowState", "maximize");
subplot(2, 1, 1);
plot(measured_current.Values, 'LineWidth', 1, 'Color', 'blue');
hold on;
plot(simulated_current.Values, 'LineWidth', 1, 'Color', 'red');
grid on;

xlabel_dict = containers.Map(["en", "ru"], ["Time, s", "Время, с"]);
ylabel_dict = containers.Map(["en", "ru"], ["Current, A", "Ток, А"]);
xlabel_translate(xlabel_dict, options.Language);
ylabel_translate(ylabel_dict, options.Language);

legend_dict = containers.Map(["en", "ru"], {["Measured", "Simulated"], ["Эксперимент", "Модель"]});
legend_translate(legend_dict, options.Language, 'Location', 'best', 'FontSize', 14);

title_dict = containers.Map(["en", "ru"], ["Motor current", "Ток двигателя"]);
title(title_dict(options.Language));

ax = gca;
ax.FontSize = 14;
ax.XLimitMethod = 'tight';
subplot(2, 1, 2);
plot(measured_velocity.Values, 'LineWidth', 1, 'Color', 'blue');
hold on;
plot(simulated_velocity.Values, 'LineWidth', 1, 'Color', 'red');
grid on;

ylabel_dict = containers.Map(["en", "ru"], ["Velocity, rad/s", "Скорость, рад/с"]);
xlabel_translate(xlabel_dict, options.Language);
ylabel_translate(ylabel_dict, options.Language);

legend_translate(legend_dict, options.Language, 'Location', 'best', 'FontSize', 14);

title_dict = containers.Map(["en", "ru"], ["Motor velocity", "Скорость двигателя"]);
title(title_dict(options.Language));

ax = gca;
ax.FontSize = 14;
ax.XLimitMethod = 'tight';
end

function resnorm = calculate_squared_norm_of_residual(simulated_dataset, measured_dataset)
arguments
    simulated_dataset Simulink.SimulationData.Dataset,
    measured_dataset Simulink.SimulationData.Dataset
end

r = sdo.requirements.SignalTracking('Method', 'Residuals', 'Normalize', 'on');

Error = [];
for ctSig=1:numElements(simulated_dataset)
    sig_name = simulated_dataset{ctSig}.Name;
    simulated_signal = getElement(simulated_dataset, sig_name);
    measured_signal = getElement(measured_dataset, sig_name);

    err = evalRequirement(r, simulated_signal.Values, measured_signal.Values);
    
    Error = [Error; err(:)];
end
resnorm = sum(Error.^2);
end