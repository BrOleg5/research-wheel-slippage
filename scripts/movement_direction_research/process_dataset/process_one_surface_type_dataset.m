%% TODO: split this into local functions
function process_one_surface_type_dataset(input_dataset_file, options)
% PROCESS_ONE_SURFACE_TYPE_DATASET Entry point

arguments
    input_dataset_file {mustBeTextScalar, mustBeNonzeroLengthText} = "output_files/datasets/for_research_movement_direction/one_surface_type/no_averaged_data.csv",
    options.KalmanFilter logical = false
end

output_dataset_file = input_dataset_file;

[robotParameters, robotKinematics] = get_robot_parameters();

% Read camera parameters
fid = fopen("config\camera_params.json");
raw = fread(fid, inf);
fclose(fid);
camera_param = jsondecode(char(raw'));
cameraParameters = struct('x', camera_param.pixel_resolution_x/1e3, ...
                          'y', camera_param.pixel_resolution_y/1e3);

check_dependency(input_dataset_file, @union_one_surface_type_dataset);

opts = detectImportOptions(input_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
opts = setvartype(opts, "surftype", "categorical");
T = readtable(input_dataset_file, opts);
% summary(T);

% Read headlines of file
fid = fopen(input_dataset_file);
column_description = fgetl(fid);
column_units = fgetl(fid);
fclose(fid);

% Remove dx, dy and dang
for var_name = ["dang", "dy", "dx"]
    idx = find(strcmp(T.Properties.VariableNames, var_name), 1);
    T = removevars(T, var_name);
    column_description = remove_item(column_description, idx);
    column_units = remove_item(column_units, idx);
end

if(options.KalmanFilter)
    % Model parameters of motor 1
    R_a = 5.95; % Ohm
    L_a = 8.9e-3; % H
    k_t = 0.0514; % N*m/A
    k_e = k_t;
    J = 1.1407e-6; % kg*m^2
    k_p = 1.6713;
    tau = 0.153; % s
    b = 2.3685e-5; % N*m*s/rad
    T_0 = 0.0081312; % N*m

    % State-space model
    % x = [i; omega; T_l];
    A = [-R_a/L_a, -(k_p + k_e)/L_a,    0; 
            k_t/J,             -b/J,  1/J;
                0,                0,    0];
    B = [k_p/L_a; 0; 0];
    C = [1, 0, 0;
         0, 1, 0];

    % Discretization
    Ts = 30e-3; % s
    Big = expm([A, B; zeros(1, 4)]*Ts);
    F = Big(1:3, 1:3);
    G = Big(1:3, 4);
    H = C;

    % Check observability of discrete system
    O = rank([H; H*F; H*(F)^2]);
    if(O == size(A, 1))
        disp("Discrete system is fully observable.");
    else
        warning("Discrete system isn't fully observable.");
    end
    
    % Process noise matrix
    Q_d = [4.83e-4,  1.02e-2,  -1.28e-5;
           1.02e-2,    26.47, -3.31e-4;
           -1.28e-5, -3.31e-4,  4.14e-7];
    if(~issymmetric(Q_d))
        error("Process noise covariance matrix Q doesn't symmetric.");
    end

    % Measurement noise matrix
    R_d = [4.80e-3,     0;
                 0,   29.57];
    if(~issymmetric(R_d))
        error("Measurement noise covariance matrix R doesn't symmetric.");
    end

    kalmanFilter = KalmanFilter(F, G, H, Q_d, R_d);
    exp = unique(T.expnum);
    for m = 1:3
        tor = zeros(height(T), 1);
        s = num2str(m);
        setvelstr = "m" + s + "setvel";
        velstr = "m" + s + "vel";
        curstr = "m" + s + "cur";
        vars = ["t", setvelstr, velstr, curstr];
        for e = exp.'
            idx = T.expnum == e;
            T1 = T(idx, vars);
            
            n = height(T1);

            filtered_vel = zeros(n, 1);
            filtered_cur = zeros(n, 1);
            estimated_tor = zeros(n, 1);

            x_hat = zeros(3, 1);
            P = [3^2,      0,       0;
                   0,  300^2,       0;
                   0,      0,   0.2^2];
            i_prev = 1;
            for i = 1:n
                if((T1.t(i) - T1.t(i_prev)) > tau)
                    i_prev = i_prev + 1;
                end
                u = T1.(setvelstr)(i_prev);
                [x_hat, P] = kalmanFilter.predict(x_hat, u, P);
                z = [T1.(curstr)(i) * sign(T1.(velstr)(i)); T1.(velstr)(i)];
                [x_hat, P] = kalmanFilter.correct(x_hat, z, P);

                filtered_cur(i) = x_hat(1);
                filtered_vel(i) = x_hat(2);
                estimated_tor(i) = x_hat(3);
            end
            T.(curstr)(idx) = abs(filtered_cur);
            T.(velstr)(idx) = filtered_vel;
            tor(idx) = estimated_tor;
        end
        torstr = "m" + s + "esttor";
        T = addvars(T, tor, 'NewVariableNames', torstr);
    end
    column_description = sprintf("%s;1st motor estimated torque;2nd motor estimated torque;3rd motor estimated torque", ...
                                 column_description);
    column_units = sprintf("%s;Nm;Nm;Nm", column_units);
end

% Calculate robot position in meter
T.xpos = T.xpos * cameraParameters.x;
idx = find(strcmp(T.Properties.VariableNames, "xpos"), 1);
column_units = replace_item(column_units, idx, "m");
T.ypos = T.ypos * cameraParameters.y;
idx = find(strcmp(T.Properties.VariableNames, "ypos"), 1);
column_units = replace_item(column_units, idx, "m");

% Calculate robot speed in local coordinate frame
[vx, vy, omega] = calc_local_robot_speed(T.t/1e3, T.xpos, T.ypos, deg2rad(T.ang));
T = addvars(T, vx, vy, omega);
column_description = sprintf("%s;X axis robot speed;Y axis robot speed;Rotational robot velocity", column_description);
column_units = sprintf("%s;m/s;m/s;rad/s", column_units);

% Calculate velocity of wheels
w1vel = T.m1vel / robotParameters.gear_ratio;
w2vel = T.m2vel / robotParameters.gear_ratio;
w3vel = T.m3vel / robotParameters.gear_ratio;
T = addvars(T, w1vel, w2vel, w3vel);
column_description = sprintf("%s;1st wheel velocity; 2nd wheel velocity;3rd wheel velocity", column_description);
column_units = sprintf("%s;rad/s;rad/s;rad/s", column_units);

% Calculate effective wheel velocities
[w1effvel, w2effvel, w3effvel] = robotKinematics.inverse(T.vx, T.vy, T.omega);
T = addvars(T, w1effvel, w2effvel, w3effvel);
column_description = sprintf("%s;1st wheel effective velocity;2nd wheel effective velocity;3rd wheel effective velocity", column_description);
column_units = sprintf("%s;rad/s;rad/s;rad/s", column_units);

% Calculate slippage of wheels
low_vel = 0.2;
w1slip = 1 - T.w1effvel ./ T.w1vel;
w1slip(abs(T.w1vel) < low_vel) = 0;
w2slip = 1 - T.w2effvel ./ T.w2vel;
w2slip(abs(T.w2vel) < low_vel) = 0;
w3slip = 1 - T.w3effvel ./ T.w3vel;
w3slip(abs(T.w3vel) < low_vel) = 0;
T = addvars(T, w1slip, w2slip, w3slip);
column_description = sprintf("%s;1st wheel slippage;2nd wheel slippage;3rd wheel slippage", column_description);
column_units = sprintf("%s;none;none;none", column_units);

% Calculate drawdown of wheel velocity
w1ddvel = T.w1vel - T.w1effvel;
w2ddvel = T.w2vel - T.w2effvel;
w3ddvel = T.w3vel - T.w3effvel;
T = addvars(T, w1ddvel, w2ddvel, w3ddvel);
column_description = sprintf("%s;1st wheel velocity drawdown;2nd wheel velocity drawdown;3rd wheel velocity drawdown", column_description);
column_units = sprintf("%s;rad/s;rad/s;rad/s", column_units);

% Calculate drawdown of robot velocity
[vx_enc, vy_enc, omega_enc] = robotKinematics.direct(T.w1vel, T.w2vel, T.w3vel);
ddvx = vx_enc - T.vx;
ddvy = vy_enc - T.vy;
ddomega = omega_enc - T.omega;
T = addvars(T, ddvx, ddvy, ddomega);
column_description = sprintf("%s;X axis robot speed drawdown;Y axis robot speed drawdown;Rotational robot velocity drawdown", column_description);
column_units = sprintf("%s;m/s;m/s;rad/s", column_units);

% Calculate drawdown of robot speed vector
v_enc = sqrt(vx_enc.^2 + vy_enc.^2);
v_eff = sqrt(T.vx.^2 + T.vy.^2);
ddspeedvec = v_enc - v_eff;
T = addvars(T, ddspeedvec);
column_description = sprintf("%s;Robot speed vector drawdown", column_description);
column_units = sprintf("%s;m/s", column_units);

% Slip linearization by taking into account the contribution of the wheel to the movement of the robot
[vx, vy, ~] = robotKinematics.direct(T.w1vel, T.w2vel, T.w3vel);
K = calc_linearization_coefficient(vx, vy, robotKinematics);
n = height(T);
w1linslip = squeeze(K(1, :, :)) .* T.w1slip;
w2linslip = squeeze(K(2, :, :)) .* T.w2slip;
w3linslip = squeeze(K(3, :, :)) .* T.w3slip;
T = addvars(T, w1linslip, w2linslip, w3linslip);
column_description = sprintf("%s;1st wheel linear slippage;2nd wheel linear slippage;3rd wheel linear slippage", column_description);
column_units = sprintf("%s;none;none;none", column_units);

% Write headlines in file
fid = fopen(output_dataset_file, "w");
fprintf(fid, "%s\n%s", column_description, column_units);
fclose(fid);
% Write table in file
writetable(T, output_dataset_file, 'Delimiter', ';', 'WriteMode', 'Append', ...
           'WriteVariableNames', true);
end
%% Local function
function [vx, vy, omega] = calc_local_robot_speed(t, X, Y, Ang)
arguments
    t {mustBeNumeric, mustBeVector, mustBeNonempty},
    X {mustBeNumeric, mustBeVector, mustBeNonempty},
    Y {mustBeNumeric, mustBeVector, mustBeNonempty},
    Ang {mustBeNumeric, mustBeVector, mustBeNonempty}
end
    assert(length(t) == length(X), "t and X vectors don't have same length.");
    assert(length(X) == length(Y), "X and Y vectors don't have same length.");
    assert(length(X) == length(Ang), "X and Ang vectors don't have same length.");
    
    dt = diff(t);
    dx = diff(X);
    dy = diff(Y);

    vx = zeros(size(X));
    vy = zeros(size(X));
    n = length(t);
    for i = 2:n
        vx(i) = (-dx(i-1) * sin(Ang(i)) - dy(i-1) * cos(Ang(i))) / dt(i-1);
        vy(i) = (-dx(i-1) * cos(Ang(i)) + dy(i-1) * sin(Ang(i))) / dt(i-1);
    end
    omega = [0; angdiff(Ang) ./ dt];
end