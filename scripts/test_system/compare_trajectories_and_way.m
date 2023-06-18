function compare_trajectories_and_way(options)
%TEST_SYSTEM Entry point

arguments
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty},
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "ru",
    options.ShowFig logical = false;
end

input_dataset_file = ...
    "output_files/datasets/for_research_movement_direction/one_surface_type/no_averaged_data.csv";

input_averaged_dataset_file = ...
    "output_files/datasets/for_research_movement_direction/one_surface_type/averaged_data.csv";

opts = detectImportOptions(input_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
T = readtable(input_dataset_file, opts);

opts = detectImportOptions(input_averaged_dataset_file, "VariableDescriptionsLine", 1, ...
                           "VariableUnitsLine", 2, "VariableNamesLine", 3, ...
                           "LeadingDelimitersRule", "error", "TrailingDelimitersRule", "ignore", ...
                           "ConsecutiveDelimitersRule", "error", "MissingRule", "error", ...
                           "ImportErrorRule", "error", "EmptyLineRule", "skip");
avrT = readtable(input_averaged_dataset_file, opts);

T = process_test_dataset(T);
avrT = process_test_dataset(avrT);

check_table_vars(T.Properties.VariableNames, ["xpos", "ypos"]);

log_folder = "output_files/logs/test_system";
if(~isfolder(log_folder))
    mkdir(log_folder);
end
log_file_name = fullfile(log_folder, "way_compare.txt");
logID = fopen(log_file_name, 'w');

head_table_format = "Exp | Surface | Speed | Direction | Camera way | Odometry way | Slip model way | Kalman + slip model way |\n";
head_separator = [repelem('-', strlength(head_table_format) - 2), '\n'];
body_table_format = "%3d | %7s | %5.1f | %9d | %10.4f | %12.4f | %14.4f | %23.4f |\n";

unique_expnum = unique(T.expnum, "sorted");
for expnum = unique_expnum.'
    idx = T.expnum == expnum;
    subT = T(idx, :);
    subT.xpos = subT.xpos - subT.xpos(1);
    subT.ypos = subT.ypos - subT.ypos(1);

    idx = avrT.expnum == expnum;
    subavrT = avrT(idx, :);
    subavrT.xpos = subavrT.xpos - subavrT.xpos(1);
    subavrT.ypos = subavrT.ypos - subavrT.ypos(1);
    
    odom_traj = odometry(subT, subT.w1vel, subT.w2vel, subT.w3vel);
    slip_model_traj = odometry_with_slip_models(subT);
    kalman_traj = odometry_kalman(subT);
    
    avr_odom_traj = odometry(subavrT, subavrT.w1vel, subavrT.w2vel, subavrT.w3vel);
    
    % Plot trajectories
    if(options.ShowFig)
        fig = figure('Name', "Test exp " + num2str(expnum));
        hold on;
        grid on;
        ax = gca;
        ax.YDir = 'reverse';
        plot(subT.xpos, subT.ypos, "Color", "black", "LineStyle", "-", "LineWidth", 1);
        plot(odom_traj.x, odom_traj.y, "Color", "black", "LineStyle", "--", "LineWidth", 1);
        plot(slip_model_traj.x, slip_model_traj.y, "Color", "black", "LineStyle", "-.", "LineWidth", 1);
        plot(kalman_traj.x, kalman_traj.y, "Color", "black", "LineStyle", ":", "LineWidth", 1);
        xlabel_dict = containers.Map(["en", "ru"], ["X, m", "X, м"]);
        ylabel_dict = containers.Map(["en", "ru"], ["Y, m", "Y, м"]);
        xlabel_translate(xlabel_dict, options.Language);
        ylabel_translate(ylabel_dict, options.Language);
        legend_dict = containers.Map(["en", "ru"], {["camera", "odometry", "odom + slip", "odom + kalman + slip"], ...
            ["камера", "одометрия", "одом + проск ", "одом + калман + проск "]});
        lgd = legend_translate(legend_dict, options.Language, 'Location', 'best');
        legend_title_dict = containers.Map(["en", "ru"], ["Trajectories", "Траектории"]);
        title(lgd, legend_title_dict(options.Language));
        % ax.XLim = [1.3, 1.7];
    end
    
    camera_traj = struct("x", subT.xpos, "y", subT.ypos);
    camera_way = calc_way(camera_traj);
    odom_way = calc_way(odom_traj);
    slip_model_way = calc_way(slip_model_traj);
    kalman_way = calc_way(kalman_traj);
    
    fprintf(logID, head_table_format);
    fprintf(logID, head_separator);
    fprintf(logID, body_table_format, subT.expnum(1), string(subT.surftype(1)), subT.speedamp(1), ...
        subT.movedir(1), camera_way, odom_way, slip_model_way, kalman_way);
    fprintf(logID, "\n");
end
fclose(logID);
end
%% Local functions
function traj = odometry(T, w1effvel, w2effvel, w3effvel)
arguments
    T table,
    w1effvel {mustBeNumeric, mustBeVector},
    w2effvel {mustBeNumeric, mustBeVector},
    w3effvel {mustBeNumeric, mustBeVector},
end

[~, robotKinematics] = get_robot_parameters();

check_table_vars(T.Properties.VariableNames, ["t", "xpos", "ypos"]);

[x_speed, y_speed, vel] = robotKinematics.direct(w1effvel, w2effvel, w3effvel);

dt = diff(T.t);
if(mean(dt) > 1)
    dt = dt / 1e3;
end

n = height(T);
traj.x = zeros(n, 1);
traj.y = zeros(n, 1);
traj.ang = zeros(n, 1);

traj.x(1) = T.xpos(1);
traj.y(1) = T.ypos(1);
traj.ang(1) = deg2rad(T.ang(1)) - pi/2;
for i = 2:n
    ang = traj.ang(i-1);

    vx = cos(ang) * x_speed(i-1) - sin(ang) * y_speed(i-1);
    traj.x(i) = traj.x(i-1) - vx * dt(i-1);

    vy = sin(ang) * x_speed(i-1) + cos(ang) * y_speed(i-1);
    traj.y(i) = traj.y(i-1) + vy * dt(i-1);

    traj.ang(i) = traj.ang(i-1) + vel(i-1) * dt(i-1);
%     traj.ang(i) = deg2rad(T.ang(i)) - pi/2;
end
end

function traj = odometry_with_slip_models(T)
arguments
    T table
end

load("output_files\models\slippage_vs_motor_linear_regression_models.mat", "lin_model");

check_table_vars(T.Properties.VariableNames, ["m1cur", "m2cur", "m3cur", ...
    "w1vel", "w2vel", "w3vel", "surftype"]);

surf_str = string(T.surftype(1));
w1slip = polyval(lin_model(1).(surf_str), abs(T.m1cur));
w2slip = polyval(lin_model(2).(surf_str), abs(T.m2cur));
w3slip = polyval(lin_model(3).(surf_str), abs(T.m3cur));

w1slip(w1slip < 0) = 0;
w2slip(w2slip < 0) = 0;
w3slip(w3slip < 0) = 0;

w1effvel = T.w1vel .* (1 - w1slip);
w2effvel = T.w2vel .* (1 - w2slip);
w3effvel = T.w3vel .* (1 - w3slip);

traj = odometry(T, w1effvel, w2effvel, w3effvel);
end

function traj = odometry_kalman(T)
arguments
    T table
end

load("output_files/model_parameters/motor_models.mat", "motor_model");

check_table_vars(T.Properties.VariableNames, ["t", "m1cur", "m2cur", "m3cur", ...
    "m1setvel", "m2setvel", "m3setvel", "m1vel", "m2vel", "m3vel"]);

P = [3^2,      0,       0;
       0,  300^2,       0;
       0,      0,   0.2^2];
x_hat = zeros(size(P, 1), 1);

for m = 1:3
    mstr = "m" + num2str(m);
    motor_data = T(:, ["t", mstr + "setvel", mstr + "cur", mstr + "vel"]);
    motor_data = renamevars(motor_data,[mstr + "setvel", mstr + "cur", mstr + "vel"], ...
        ["setvel", "cur", "vel"]);

    x_hat_history = kalman_filter_process(x_hat, P, motor_data, motor_model(m));
    
    wstr = "w" + num2str(m);
    T.(mstr + "cur") = squeeze(x_hat_history(1, 1, :));
    T.(mstr + "vel") = squeeze(x_hat_history(2, 1, :));
    T.(wstr + "vel") = T.(mstr + "vel") / 16;
end

traj = odometry_with_slip_models(T);
end

function way = calc_way(traj)
arguments
    traj struct
end

assert(all(isfield(traj, ["x", "y"])), "Struct must have this fields");

dx = diff(traj.x);
dy = diff(traj.y);
way = sum(sqrt(dx.^2 + dy.^2));
end