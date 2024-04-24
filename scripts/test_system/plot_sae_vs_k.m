function plot_sae_vs_k(options)
arguments
    options.ExportGraphs (1, 1) {mustBeNumericOrLogical} = false,
    options.ExportGraphExtensions {mustBeText, mustBeNonzeroLengthText, mustBeNonempty, ...
        mustBeMember(options.ExportGraphExtensions, ["jpg", "jpeg", "png", "tif", "tiff", "gif", ...
        "eps", "emf", "pdf"])} = ["emf", "pdf"]
    options.ExportGraphFolder {mustBeText, mustBeNonempty} = "output_files/graphs/test_system",
    options.Language {mustBeTextScalar, mustBeMember(options.Language, ["en", "ru"])} = "en"
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

xlabel_dict = containers.Map(["en", "ru"], ["Time, s", "Время, с"]);
ylabel_dict = containers.Map(["en", "ru"], ["Velocity, rad/s", "Угловая скорсть, рад/с"]);
legend_dict = containers.Map(["en", "ru"], {["slippage model", "camera"], ...
    ["модель проскальзывания", "камера"]});

speedamp = 0.1;
idx = T.speedamp == speedamp;
T = T(idx, :);
idx = avrT.speedamp == speedamp;
avrT = avrT(idx, :);

unique_expnum = unique(T.expnum, "sorted");
n = length(unique_expnum);

varTypes = repelem("double", 1, 12);
result = table('Size', [n, 12], 'VariableTypes', varTypes, ...
    'VariableNames', ["k1", "k2", "k3", "w1slip", "w2slip", "w3slip", "w1linslip", "w2linslip", ...
    "w3linslip", "mse1", "mse2", "mse3"]);
[~, robotKinematics] = get_robot_parameters();
for i = 1:n
    expnum = unique_expnum(i);
    idx = T.expnum == expnum;
    subT = T(idx, :);

    idx = avrT.expnum == expnum;
    subavrT = avrT(idx, :);

    [w1vel, w2vel, w3vel, ~, ~] = kalman_f(subT);
    kalman = struct("w1vel", w1vel, "w2vel", w2vel, "w3vel", w3vel);

    % Interpolate and extrapolate wheel effective velocity in averaged dataset
%     interp_vel = struct("w1vel", zeros(size(w1vel)), "w2vel", zeros(size(w2vel)), ...
%         "w3vel", zeros(size(w3vel)));
%     for w = 1:3
%         wstr = "w" + num2str(w) + "vel";
%         weffstr = "w" + num2str(w) + "effvel";
%         interp_vel.(wstr) = interp1(subavrT.t.*1e3, subavrT.(weffstr), subT.t, "linear", "extrap");
%     end
% 
%     result.mse1(i) = mean((kalman.w1vel - interp_vel.w1vel).^2);
%     result.mse2(i) = mean((kalman.w2vel - interp_vel.w2vel).^2);
%     result.mse3(i) = mean((kalman.w3vel - interp_vel.w3vel).^2);

    result.mse1(i) = mean((kalman.w1vel - subT.w1effvel).^2);
    result.mse2(i) = mean((kalman.w2vel - subT.w2effvel).^2);
    result.mse3(i) = mean((kalman.w3vel - subT.w3effvel).^2);
    
    assert(all(~diff(subT.movedir)), "All movement direction must be equal in one experiment");
    assert(all(~diff(subavrT.movedir)), "All movement direction must be equal in one experiment");
    assert(subT.movedir(1) == subavrT.movedir(1), "Movement direction in subT and subavrT must be equal");
    K = calc_linearization_coefficient(deg2rad(subT.movedir(1)), robotKinematics);
    
    result.k1(i) = K(1);
    result.k2(i) = K(2);
    result.k3(i) = K(3);

    result.w1slip(i) = mean(subT.w1slip);
    result.w2slip(i) = mean(subT.w2slip);
    result.w3slip(i) = mean(subT.w3slip);

    result.w1linslip(i) = mean(subT.w1linslip);
    result.w2linslip(i) = mean(subT.w2linslip);
    result.w3linslip(i) = mean(subT.w3linslip);
end


for m = 1:3
    fig = figure("Name", "Test " + num2str(m), "WindowState", "maximized");
    tiledlayout(3, 1, "TileSpacing", "tight", "Padding", "compact");
    nexttile;
    hold on;
    grid on;
    msestrstr = "mse" + num2str(m);
    kstr = "k" + num2str(m);
    plot(result.(kstr), result.(msestrstr), "LineStyle", "none", "Marker", ".", "MarkerSize", 8, ...
        "MarkerEdgeColor", "black");

    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    if(m == 2)
        legend_translate(legend_dict, options.Language, 'Location', 'best');
    end
    
    ax = gca;
    ax.FontSize = 12;
    ax.XLimitMethod = "tight";
    ax.YLimitMethod = "padded";


    nexttile;
    hold on;
    grid on;
    slipstr = "w" + num2str(m) + "slip";
    plot(result.(slipstr), result.(msestrstr), "LineStyle", "none", "Marker", ".", "MarkerSize", 8, ...
        "MarkerEdgeColor", "black");

    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    if(m == 2)
        legend_translate(legend_dict, options.Language, 'Location', 'best');
    end
    
    ax = gca;
    ax.FontSize = 12;
    ax.XLimitMethod = "tight";
    ax.YLimitMethod = "padded";


    nexttile;
    hold on;
    grid on;
    linslipstr = "w" + num2str(m) + "linslip";
    plot(result.(linslipstr), result.(msestrstr), "LineStyle", "none", "Marker", ".", "MarkerSize", 8, ...
        "MarkerEdgeColor", "black");

    xlabel_translate(xlabel_dict, options.Language);
    ylabel_translate(ylabel_dict, options.Language);
    
    ax = gca;
    ax.FontSize = 12;
    ax.XLimitMethod = "tight";
    ax.YLimitMethod = "padded";
    if(options.ExportGraphs)
        output_dir = fullfile(options.ExportGraphFolder, subT.surftype(1));
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, num2str(subT.speedamp(1)) + " " + num2str(subT.movedir(1)));
        export_graphs(fig, file_name, options.ExportGraphExtensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end
%% Local functions
function [w1vel, w2vel, w3vel] = slip_models(T)
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

w1vel = T.w1vel .* (1 - w1slip);
w2vel = T.w2vel .* (1 - w2slip);
w3vel = T.w3vel .* (1 - w3slip);
end

function [w1vel, w2vel, w3vel, kalman_cur, kalman_vel] = kalman_f(T)
arguments
    T table
end

load("output_files/model_parameters/motor_models.mat", "motor_model");

check_table_vars(T.Properties.VariableNames, ["t", "m1cur", "m2cur", "m3cur", ...
    "surftype", "m1setvel", "m2setvel", "m3setvel", "m1vel", "m2vel", "m3vel", "movedir"]);

P_0 = [3^2,      0,       0;
       0,  300^2,       0;
       0,      0,   0.2^2];
x_hat_0 = zeros(size(P_0, 1), 1);

kalman_cur = {};
kalman_vel = {};

filtered_T = table;
filtered_T.surftype = T.surftype;
filtered_T.movedir = T.movedir;
for m = 1:3
    mstr = "m" + num2str(m);
    wstr = "w" + num2str(m);
    motor_data = T(:, ["t", mstr + "setvel", mstr + "cur", mstr + "vel"]);
    motor_data = renamevars(motor_data, [mstr + "setvel", mstr + "cur", mstr + "vel"], ...
        ["setvel", "cur", "vel"]);

    x_hat = kalman_filter_process(x_hat_0, P_0, motor_data, motor_model(m));

    kalman_cur{m} = squeeze(x_hat(1, 1, :));
    kalman_vel{m} = squeeze(x_hat(2, 1, :)) / 16;

    filtered_T = addvars(filtered_T, kalman_cur{m}, kalman_vel{m}, ...
        'NewVariableNames', [mstr + "cur", wstr + "vel"]);
end

[w1vel, w2vel, w3vel] = slip_models(filtered_T);
% [w1vel, w2vel, w3vel] = slip_nonlinear_models(fisltered_T);
end

function [w1vel, w2vel, w3vel] = slip_nonlinear_models(T)
arguments
    T table
end

[~, robotKinematics] = get_robot_parameters();

load("output_files\models\slippage_vs_motor_linear_regression_models.mat", "lin_model");

check_table_vars(T.Properties.VariableNames, ["m1cur", "m2cur", "m3cur", ...
    "w1vel", "w2vel", "w3vel", "surftype", "movedir"]);

surf_str = string(T.surftype(1));
w1slip = polyval(lin_model(1).(surf_str), abs(T.m1cur));
w2slip = polyval(lin_model(2).(surf_str), abs(T.m2cur));
w3slip = polyval(lin_model(3).(surf_str), abs(T.m3cur));

w1slip(w1slip < 0) = 0;
w2slip(w2slip < 0) = 0;
w3slip(w3slip < 0) = 0;

K = calc_linearization_coefficient(T.movedir, robotKinematics);

w1slip = w1slip ./ squeeze(K(1, :, :));
w2slip = w2slip ./ squeeze(K(2, :, :));
w3slip = w3slip ./ squeeze(K(3, :, :));

w1vel = T.w1vel .* (1 - w1slip);
w2vel = T.w2vel .* (1 - w2slip);
w3vel = T.w3vel .* (1 - w3slip);
end