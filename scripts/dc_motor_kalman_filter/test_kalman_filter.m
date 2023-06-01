function test_kalman_filter()

export_graph = true;
export_graph_folder = "output_files/graphs/dc_motor_kalman_filter";
export_graph_extensions = ["emf", "pdf"];

motor_models_path = "output_files/model_parameters/motor_models.mat";
load(motor_models_path, "motor_model");

% Load dataset
dataset = readtable("datasets\dc_motor_model_identification\data_for_identification_without_integral.csv", ...
                       "NumHeaderLines", 2);
% dataset = readtable("datasets\dc_motor_model_identification\validation_data.csv", ...
%                        "NumHeaderLines", 2);
% dataset = readtable("datasets\dc_motor_model_identification\green_gray_validation_data.csv", ...
%                        "NumHeaderLines", 2);
% dataset = readtable("datasets\dc_motor_model_identification\table_green_validation_data.csv", ...
%                        "NumHeaderLines", 2);
dataset = process_identification_dataset(dataset);

P = [3^2,      0,       0;
       0,  300^2,       0;
       0,      0,   0.2^2];
x_hat = zeros(size(P, 1), 1);

for motor_num = 1:3
    mstr = "m" + num2str(motor_num);
    motor_data = dataset(:, ["t", mstr + "setvel", mstr + "cur", mstr + "vel"]);
    motor_data = renamevars(motor_data,[mstr + "setvel", mstr + "cur", mstr + "vel"], ...
        ["setvel", "cur", "vel"]);

    [x_hat_history, P_history] = simulate(x_hat, P, motor_data, motor_model(motor_num));

    ci = calc_confidence_interval(x_hat_history, P_history);
    ci = [squeeze(ci(:, 1, :)), squeeze(ci(:, 2, end:-1:1))];

    fig = plot_compare(x_hat_history, P_history, motor_data, ci);

    if(export_graph)
        output_dir = fullfile(export_graph_folder, "motor" + num2str(motor_num));
        if(~isfolder(output_dir))
            mkdir(output_dir);
        end
        file_name = fullfile(output_dir, "train_data");
        export_graphs(fig, file_name, export_graph_extensions, 'ContentType', 'vector', ...
                      'BackgroundColor', 'white', 'Colorspace', 'rgb');
        close(fig);
    end
end
end
%% Local functions
function [x_hat_history, P_history] = simulate(x_hat, P, dataset, model, options)
arguments
    x_hat (3, 1) {mustBeNumeric, mustBeReal},
    P (3, 3) {mustBeNumeric, mustBeReal},
    dataset table,
    model struct,
    options.ReferenceDelay = 140e-3
end

check_table_vars(dataset.Properties.VariableNames, ["t", "setvel", "cur", "vel"]);

i_prev = 1;
x_hat_history = x_hat;
P_history = P;
for i = 1:height(dataset)
    if((dataset.t(i) - dataset.t(i_prev)) > options.ReferenceDelay)
        i_prev = i_prev + 1;
    end
    [x_hat, P] = predict(x_hat, dataset.setvel(i_prev), P, model);
    z = [dataset.cur(i); dataset.vel(i)];
    [x_hat, P] = correct(x_hat, z, P, model);
    x_hat_history(:, :, i) = x_hat;
    P_history(:, :, i) = P;
end
end

function fig = plot_compare(x_hat_history, P_history, dataset, ci)
arguments
    x_hat_history (3, 1, :) {mustBeNumeric, mustBeReal},
    P_history (3, 3, :) {mustBeNumeric, mustBeReal},
    dataset table,
    ci {mustBeNumeric, mustBeReal}
end

n_x = size(x_hat_history, 1);

check_table_vars(dataset.Properties.VariableNames, ["t", "setvel", "cur", "vel"]);

fig = figure("Name", "Kalman filter", "WindowState", "maximized");
y_name = ["Ток, А", "Скорость, рад/с", "Момент, Нм"];
t_name = ["cur", "vel"];
uns_name = ["Неопределённость тока", "Неопределённость скорости", "Неопределённость момента"];
tiledlayout(n_x, 2, "TileSpacing", "tight", "Padding", "compact");
for i = 1:n_x
    nexttile;
    fill([dataset.t; dataset.t(end:-1:1)], ci(i, :), "green", "FaceAlpha", 0.2, "EdgeColor", "none");
    grid on;
    hold on;
    if(i ~= 3)
        plot(dataset.t, dataset.(t_name(i)), 'Color', 'blue', 'LineWidth', 0.75);
    end
    plot(dataset.t, squeeze(x_hat_history(i, :, :)), 'Color', 'red', 'LineWidth', 0.75);
    xlabel("Время, с");
    ylabel(y_name(i));
    if(i == 2)
        legend(["Доверительный интервал 95%", "Измерения", "Отфильтрованное"], ...
            "Location", "best", "FontSize", 11);
    end
    ax = gca;
    ax.FontSize = 11;
    nexttile;
    plot(dataset.t, squeeze(P_history(i, i, :)), 'Color', 'black', 'LineWidth', 0.75);
    grid on;
    hold on;
    xlabel("Время, с");
    ylabel(uns_name(i));
    ax = gca;
    ax.FontSize = 11;
end
end

% Kalman filter function
function [x_hat, P] = predict(x_hat, u, P, model)
arguments
    x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
    u (:, 1) {mustBeNumeric, mustBeNonempty},
    P (:, :) {mustBeNumeric, mustBeNonempty},
    model struct {mustBeNonempty}
end
    assert(size(P, 1) == size(P, 2), "P matrix must be square.");
    assert(size(model.F, 1) == size(model.F, 2), "F matrix must be square.");
    assert(size(model.Q, 1) == size(model.Q, 2), "Q matrix must be square.");
    assert(length(x_hat) == size(P, 1), "Number of state x_hat must be equal size of P matrix.");
    assert(length(x_hat) == size(model.Q, 1), "Number of state x_hat must be equal size of Q matrix.");
    assert(length(x_hat) == size(model.F, 1), "Number of state x_hat must be equal size of F matrix.");
    assert(length(x_hat) == size(model.G, 1), "Number of state x_hat  must be equal row number of G matrix.");
    assert(length(u) == size(model.G, 2), "Number of control variable u must be equal column number of G matrix.");
    
    x_hat = model.F * x_hat + model.G * u;
    P = model.F * P * model.F.' + model.Q;
end

function [x_hat, P, K] = correct(x_hat, z, P, model)
arguments
    x_hat (:, 1) {mustBeNumeric, mustBeNonempty},
    z (:, 1) {mustBeNumeric, mustBeNonempty},
    P (:, :) {mustBeNumeric, mustBeNonempty},
    model struct {mustBeNonempty}
end
    assert(size(P, 1) == size(P, 2), "P matrix must be square.");
    assert(size(model.R, 1) == size(model.R, 2), "R matrix must be square.");
    assert(length(x_hat) == size(P, 1), "Number of state x_hat must be equal size of P matrix.");
    assert(length(x_hat) == size(model.H, 2), "Number of state x_hat must be equal column number of H matrix.");
    assert(length(z) == size(model.R, 1), "Number of measurements z must be equal size of R matrix.");
    assert(length(z) == size(model.H, 1), "Number of measurements z must be equal row number of H matrix.");

    K = P * model.H.' / (model.H * P * model.H.' + model.R);
    x_hat = x_hat + K * (z - model.H * x_hat);
    I = eye(length(x_hat));
    Z = I - K * model.H;
%     P = Z * P * Z.' + K * model.R * K.';
    P = Z * P;
end

% Confidence interval
function ci = calc_confidence_interval(x_hat, P)
arguments
    x_hat (:, :, :) {mustBeNumeric, mustBeNonempty}
    P (:, :, :) {mustBeNumeric, mustBeNonempty}
end
    for i = 1:size(x_hat, 3)
        for j = 1:size(x_hat, 1)
            if(P(j, j, i) ~= 0)
                ci(j, :, i) = norminv([0.025, 0.975], x_hat(j, 1, i), sqrt(P(j, j, i)));
            else
                ci(j, :, i) = x_hat(j, 1, i) * [1, 1];
            end
        end
    end
end