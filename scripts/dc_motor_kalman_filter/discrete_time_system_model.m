function discrete_time_system_model(options)
arguments
    options.Ts (1, 1) {mustBeNumeric, mustBeReal} = 30e-3; % s
end

model_parameters_path = "output_files/model_parameters";

% Load dataset
train_data = readtable("datasets\dc_motor_model_identification\data_for_identification_without_integral.csv", ...
                       "NumHeaderLines", 2);
train_data = process_identification_dataset(train_data);
train_data = sync_reference_and_measurements(train_data);

motor_model = struct('F', [], 'G', [], 'H', [], 'Q', [], 'R', []);
for motor_num = 1:3
    load(fullfile(model_parameters_path, "motor" + num2str(motor_num) + ".mat"), "vOpt");
    motor_param = param2struct(vOpt);
    
    % Continuous state space model
    % x = [i; omega; T_l];
    A = [-motor_param.R_a/motor_param.L_a, -(motor_param.k_p + motor_param.k_e)/motor_param.L_a,               0; 
            motor_param.k_t/motor_param.J,                         -motor_param.b/motor_param.J, 1/motor_param.J;
                                        0,                                                     0,              0];
    B = [motor_param.k_p/motor_param.L_a; 0; 0];
    C = [1, 0, 0;
         0, 1, 0];

    % Discretization
    Big = expm([A, B; zeros(1, 4)]*options.Ts);
    F = Big(1:3, 1:3);
    G = Big(1:3, 4);
    % Observation matrix
    H = C;

    % Check observability of discrete system
    O = rank([H; H*F; H*(F)^2]);
    if(O == size(A, 1))
        disp("Discrete system is fully observable.");
    else
        warning("Discrete system isn't fully observable.");
    end

    % Constructing the process noise matrix Q and measurement noise matrix R
    % Autocovariance Least-Squares
    mstr = "m" + num2str(motor_num);
    u = train_data.(mstr + "setvel").';
    z = [train_data.(mstr + "cur").'; train_data.(mstr + "vel").'];

    [data, model, estimator] = setup_als(F, G, H, u, z);
    N = 30;
    [Qest_cell, Rest_cell] = als_sdp_mrQ(data, N, model, estimator);
    Q_d = Qest_cell{1};
    R_d = Rest_cell{1};
    % Manual correction
    Q_d(1, 1) = Q_d(1, 1) / 10;
    Q_d(1, 2) = Q_d(1, 2) / 10;
    Q_d(2, 1) = Q_d(2, 1) / 10;
    Q_d(2, 3) = Q_d(2, 3) / 10;
    Q_d(3, 2) = Q_d(3, 2) / 10;

    R_d(1, 1) = R_d(1, 1) * 100;
    R_d(1, 2) = 0;
    R_d(2, 1) = 0;

    if(~issymmetric(Q_d))
        error("Process noise covariance matrix Q doesn't symmetric.");
    end

    if(~issymmetric(R_d))
        error("Measurement noise covariance matrix R doesn't symmetric.");
    end

    motor_model(motor_num) = struct('F', F, 'G', G, 'H', H, 'Q', Q_d, 'R', R_d);
end
model_parameters_path = "output_files/model_parameters";
if(~isfolder(model_parameters_path))
    mkdir(model_parameters_path);
end
save(fullfile(model_parameters_path, "motor_models.mat"), "motor_model");
end
%% Local functions
function out_struct = param2struct(params)
arguments
    params param.Continuous
end

out_struct = struct;
for i = 1:length(params)
    out_struct.(params(i).Name) = params(i).Value;
end
end

function [data, model, estimator] = setup_als(F, G, H, u, z)
arguments
    F (3, 3) {mustBeNumeric, mustBeReal},
    G (3, 1) {mustBeNumeric, mustBeReal},
    H (2, 3) {mustBeNumeric, mustBeReal},
    u (1, :) {mustBeNumeric, mustBeReal},
    z (2, :) {mustBeNumeric, mustBeReal}
end
% initial guesses
% Encoder measurement noise
% link: https://www.researchgate.net/figure/Diagram-of-speed-measurement-noise-caused-by-encoder-pulse-counting-rotor-speed-50-rpm_fig1_347615973
encoder_resolution = 500 * 2 * 2; % ticks per revolition (including two channel and rising and falling)
f_mot = 900; % Hz
T_mot_s = 1 / f_mot;
omega_noise = 2 * pi / (encoder_resolution * T_mot_s); % rad/s
omega_noise_sigma = omega_noise / 6;
omega_variance = omega_noise_sigma.^2;

current_variance = (1.5/156/6).^2;

G_hat = eye(3);
Qw_hat = [1e-7,      0,       0;
             0,   2e-1,       0;
             0,      0,    1e-9];
Rv_hat= [current_variance,                0;
                        0,   omega_variance];

L = dlqe(F, G_hat, H, Qw_hat, Rv_hat);

model.A = F;
model.B = G;
model.C = H;
model.G = G_hat;
model.xhat0 = [z(:, 1); 0];

data.datapts = length(u);
data.yk = z;
data.uk = u;

estimator.L = L;
estimator.Q = Qw_hat;
estimator.R = Rv_hat;
end